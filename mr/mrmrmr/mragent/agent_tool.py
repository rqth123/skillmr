import json
from Bio import Entrez
import requests
from bs4 import BeautifulSoup
import pandas as pd
#from rpy2.robjects import r
import time
import subprocess
import os
import json
import re
import time
import urllib


def timer(func):
    def func_wrapper(*args, **kwargs):
        from time import time
        time_start = time()
        result = func(*args, **kwargs)
        time_end = time()
        time_spend = time_end - time_start
        print('%s cost time: %.3f s' % (func.__name__, time_spend))
        return result

    return func_wrapper


def retry_with_backoff(func, retries=5, base_delay=1, max_delay=32):
    """
    Exponential backoff retry for API calls.
    Args:
        func: function to call (must accept no args)
        retries: maximum number of retries
        base_delay: base delay in seconds
        max_delay: maximum delay in seconds
    Returns:
        Result of the function call
    Raises:
        Exception if all retries fail
    """
    import time
    last_exception = None
    for attempt in range(retries):
        try:
            return func()
        except Exception as e:
            last_exception = e
            delay = min(base_delay * (2 ** attempt), max_delay)
            print(f"[retry] Attempt {attempt + 1}/{retries} failed: {e}. Retrying in {delay} seconds...")
            time.sleep(delay)
    raise last_exception


# 爬虫爬取PubMed数据
def pubmed_crawler(keyword, num_records, sort_order, json_str=True):
    Entrez.email = 'lyjjj@tmu.edu.cn'  # 请替换为你的邮箱

    def search_pubmed(keyword, num_records, sort_order):
        # 使用内嵌函数来适配无参的 retry_with_backoff
        def _do_search():
            handle = Entrez.esearch(db='pubmed',
                                    sort=sort_order,
                                    retmax=str(num_records),
                                    retmode='xml',
                                    term=keyword)
            results = Entrez.read(handle)
            return results
        
        # 增加重试机制，避免偶尔的网络超时/中断
        return retry_with_backoff(_do_search, retries=5)

    def fetch_details(id_list):
        ids = ','.join(id_list)
        
        # 使用内嵌函数包装 efetch 请求
        def _do_fetch():
            handle = Entrez.efetch(db='pubmed',
                                   retmode='xml',
                                   id=ids)
            results = Entrez.read(handle, validate=False)
            return results
            
        # 增加重试机制，防止 IncompleteRead 异常
        return retry_with_backoff(_do_fetch, retries=5)

    def get_paper_details(paper):
        try:
            title = paper['MedlineCitation']['Article']['ArticleTitle']
        except KeyError:
            title = 'NULL'

        try:
            abstract = paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
        except KeyError:
            abstract = 'NULL'

        return title, abstract

    def search_and_print_papers(keyword, num_records, sort_order):
        results = search_pubmed(keyword, num_records, sort_order)
        # print(results)
        id_list = results['IdList']
        if len(id_list) == 0:
            return json.dumps([{
                'index': 0,
                'title': 'No paper found',
                'abstract': 'No paper found'
            }])
        papers = fetch_details(id_list)
        # print(papers)
        paper_details = []
        for i, paper in enumerate(papers['PubmedArticle'], start=1):
            title, abstract = get_paper_details(paper)
            paper_details.append({
                'index': i,
                'title': title,
                'abstract': abstract
            })
        if json_str:
            return json.dumps(paper_details)
        else:
            return paper_details

    return search_and_print_papers(keyword, num_records, sort_order)
def get_gwas_id(keyword):
    def _do_request():
        url = "https://api.opengwas.io/api/gwasinfo"
        
        # 同样在这里加上你的 Token
        headers = {
            'Authorization': 'Bearer eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIzMDYwMTAwMThAcXEuY29tIiwiaWF0IjoxNzc2MTA0MDI5LCJleHAiOjE3NzczMTM2Mjl9.YMZOORxrRl28mq64xMtNDQbVbvJNBLzbwg2kHaqR3TgjKJg4alJHGsrXsOZTIPpfIB0iECzrpAOHA7u6sNY2m1VuE1Fle4KG6W2O0J2XYW7QKJ-PF7AjXw6CjcLK3uKrvunz50bmHa4dPsJwaGOvJIvs0fOULxRxN6oBa7RJb0u7LYk0heqAAed-0UDiklOGLRFJDqptw1hTJwgNTDXJ5PN4kRWLnC0eyblzNz3icKHgNuqgceQunYkopE9uTft5oUl52SguKuKkpkjnWE3Qt54OdLl-WUQivpIik0UHfeGYlgHrwSSI2HlnQbfdJ3DRVpy5YqivQtkB6Dv5CIISJw'
        }
        
        # 在 requests.get 中带上 headers
        response = requests.get(url, headers=headers, timeout=15)
        response.raise_for_status()
        data = response.json()
        
        results = []
        for d in data:
            if isinstance(d, dict) and 'trait' in d and d['trait']:
                if keyword.lower() in d['trait'].lower():
                    dataset_info = {
                        "gwas_id": d.get("id"),
                        "trait": d.get("trait"),
                        "year": d.get("year"),
                        "consortium": d.get("consortium"),
                        "Sample size": d.get("sample_size"),
                        "Number of SNPs": d.get("nsnp")
                    }
                    results.append(dataset_info)
                    
        # 数据去重
        unique_results = [dict(t) for t in {tuple(d.items()) for d in results}]
        
        # 限制返回数量，防止大模型 token 爆炸
        if len(unique_results) > 30:
            unique_results = unique_results[:30]
            
        print(f"Found {len(unique_results)} records online for '{keyword}'.")
        return unique_results

    try:
        return retry_with_backoff(_do_request, retries=3)
    except Exception as e:
        print(f"Error fetching GWAS IDs from OpenGWAS API v4: {e}")
        return []
# 获取PubMed数据库中某个文章的详细信息
def get_paper_details(paper_title):
    Entrez.email = 'lyjjj@tmu.edu.cn'  # 请替换为你的邮箱

    def search_pubmed(keyword, num_records, sort_order):
        handle = Entrez.esearch(db='pubmed',
                                sort=sort_order,
                                retmax=str(num_records),
                                retmode='xml',
                                term=keyword)
        results = Entrez.read(handle)
        return results

    def fetch_details(id_list):
        ids = ','.join(id_list)
        handle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                               id=ids)
        results = Entrez.read(handle)
        return results

    def search_and_print_papers(keyword):
        results = search_pubmed(keyword, 1, 'relevance')
        # print(results)
        id_list = results['IdList']
        papers = fetch_details(id_list)
        print(papers)
        paper = papers['PubmedArticle'][0]

        try:
            title = paper['MedlineCitation']['Article']['ArticleTitle']
        except KeyError:
            title = 'NULL'

        try:
            abstract = paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
        except KeyError:
            abstract = 'NULL'

        try:
            PMID = paper['MedlineCitation']['PMID']
        except KeyError:
            PMID = 'NULL'

        try:
            DOI = paper['MedlineCitation']['Article']['ELocationID'][0]
        except KeyError:
            DOI = 'NULL'

        try:
            Journal = paper['MedlineCitation']['Article']['Journal']['Title']
        except KeyError:
            Journal = 'NULL'

        try:
            Author = paper['MedlineCitation']['Article']['AuthorList'][0]['LastName']
        except KeyError:
            Author = 'NULL'

        try:
            Year = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
        except KeyError:
            Year = 'NULL'

        print(title)
        print(PMID)
        print(DOI)
        print(Journal)
        print(Author)
        print(Year)
        print(abstract)
        return [title, PMID, DOI, Journal, Author, Year, abstract]

    return search_and_print_papers(paper_title)


# 检测关键词是否在OpenGWAS数据库中有对应的SNP
def check_keyword_in_opengwas(keyword):
    # OpenGWAS网页地址
    def _do_request():
        url = f"https://gwas.mrcieu.ac.uk/datasets/?trait__icontains={keyword}"
        # 发送GET请求
        response = requests.get(url, timeout=30)
        # 检查响应状态码
        if response.status_code != 200:
            response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        # 检查是否存在指定的文本
        if "Filtered to 0 records" in soup.text:
            print(f"No SNP related to '{keyword}' was found in the OpenGWAS database.")
            return False
        else:
            print(f"SNP(s) related to '{keyword}' exist in the OpenGWAS database.")
            return True

    try:
        return retry_with_backoff(_do_request, retries=3)
    except Exception as e:
        print(f"An error occurred while querying the OpenGWAS database: {e}")
        return None



def get_paper_details_pmc(paper_title):
    Entrez.email = 'lyjjj@tmu.edu.cn'  # 请替换为你的邮箱

    def search_pubmed(keyword, num_records, sort_order):
        handle = Entrez.esearch(db='pubmed',
                                sort=sort_order,
                                retmax=str(num_records),
                                retmode='xml',
                                term=keyword)
        results = Entrez.read(handle)
        return results

    def get_pcm_full_text(pmc_id):
        # 直接访问接口网站获取json全文
        # https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/[ID]/unicode
        url = f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{pmc_id}/unicode"
        print(url)
        response = urllib.request.urlopen(url)
        full_text = response.read().decode('utf-8')
        return full_text

    def search_and_print_papers(keyword):
        results = search_pubmed(keyword, 1, 'relevance')
        id_list = results['IdList']
        if len(id_list) == 0:
            print("No papers found for the given query.")
            return None
        else:
            paper_id = id_list[0]
            full_text = get_pcm_full_text(paper_id)
            if 'No record can be found for the input' in full_text:
                print("No full text found for the given query.")
                return None
            else:
                return full_text

    return search_and_print_papers(paper_title)


@timer
def MRtool(Exposure_id, Outcome_id, path, gwas_token):
    """
    Run standard MR analysis using standalone R script.
    Uses subprocess for clean separation and proper error handling.
    """
    import subprocess
    import sys

    # Get the path to the R script relative to this file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path = os.path.join(script_dir, 'r_scripts', 'run_mr.R')

    # Build the command
    cmd = [
        'Rscript',
        '--vanilla',
        r_script_path,
        Exposure_id,
        Outcome_id,
        path,
        str(gwas_token) if gwas_token else ""
    ]

    print(f"[MRtool] Running command: {' '.join(cmd)}")

    # Wait before executing to avoid rate limiting
    time.sleep(5)

    # Run the process and capture output
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Print output for debugging
    if result.stdout:
        print("[MRtool] R stdout:\n", result.stdout)
    if result.stderr:
        print("[MRtool] R stderr:\n", result.stderr)

    if result.returncode != 0:
        raise RuntimeError(f"MR analysis failed with exit code {result.returncode}. Check R output above.")

    print("[MRtool] MR analysis completed successfully")


@timer
def MRtool_MOE(Exposure_id, Outcome_id, path, gwas_token):
    """
    Run MR-Mixture of Experts analysis using standalone R script.
    Uses subprocess for clean separation and proper error handling.
    """
    import subprocess
    import sys

    # Get the path to the R script relative to this file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path = os.path.join(script_dir, 'r_scripts', 'run_mr_moe.R')

    # Build the command
    cmd = [
        'Rscript',
        '--vanilla',
        r_script_path,
        Exposure_id,
        Outcome_id,
        path,
        str(gwas_token) if gwas_token else ""
    ]

    print(f"[MRtool_MOE] Running command: {' '.join(cmd)}")

    # Wait before executing to avoid rate limiting
    time.sleep(5)

    # Run the process and capture output
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Print output for debugging
    if result.stdout:
        print("[MRtool_MOE] R stdout:\n", result.stdout)
    if result.stderr:
        print("[MRtool_MOE] R stderr:\n", result.stderr)

    if result.returncode != 0:
        raise RuntimeError(f"MR-MOE analysis failed with exit code {result.returncode}. Check R output above.")

    print("[MRtool_MOE] MR-MOE analysis completed successfully")


@timer
def MRtool_MRlap(Exposure_id, Outcome_id, path, N_exposure, N_outcome):
        r_script = """
# 安装并加载必要的包
if (!requireNamespace("httr", quietly = TRUE)) {
    install.packages("httr")
}
if (!requireNamespace("vcfR", quietly = TRUE)) {
    install.packages("vcfR")
}
if (!requireNamespace("MRlap", quietly = TRUE)) {
    install.packages("MRlap")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite")
}

library(httr)
library(vcfR)
library(MRlap)
library(jsonlite)

# 定义下载函数
download_vcf <- function(base_url, file_name, dest_dir = ".", min_size = 1 * 1024 * 1024) {
    # 提取文件标识符 (如 ukb-b-10807)
    file_id <- gsub("(.*)\\.vcf\\.gz", "\\1", file_name)

    # 构建完整的文件 URL (如 ukb-b-10807/ukb-b-10807.vcf.gz)
    file_url <- paste0(base_url, file_id, "/", file_name)
    print(paste("Attempting to download from:", file_url))

    # 构建文件的本地保存路径
    dest_file <- file.path(dest_dir, file_name)
    print(paste("Saving to:", dest_file))

    # 检查文件是否已经存在并且大小是否合适
    if (file.exists(dest_file)) {
        file_info <- file.info(dest_file)
        file_size <- file_info$size

        # 如果文件存在且大小大于 1MB，则不重新下载
        if (file_size >= min_size) {
            message(paste("File already exists and is larger than", round(min_size / (1024 * 1024), 2), "MB. Skipping download."))
            return(TRUE)  # 文件已经存在且大小合适，跳过下载
        } else {
            message(paste("File exists but is too small (", round(file_size / (1024 * 1024), 2), "MB). Redownloading..."))
        }
    }

    # 尝试直接下载文件，并捕获可能的错误
    tryCatch({
        # 使用 curl 来下载并显示进度条
        download.file(file_url, destfile = dest_file, method = "curl", mode = "wb")

        # 下载完成后检查文件大小
        file_info <- file.info(dest_file)
        file_size <- file_info$size

        # 如果文件小于 1MB，视为下载失败
        if (file_size < min_size) {
            message(paste("File size is too small (", round(file_size / (1024 * 1024), 2), "MB). Deleting file:", dest_file))
            file.remove(dest_file)  # 删除文件
            return(FALSE)  # 返回 FALSE 表示下载失败
        } else {
            message(paste("Download complete:", dest_file, "File size:", round(file_size / (1024 * 1024), 2), "MB"))
            return(TRUE)  # 返回 TRUE 表示下载成功
        }

    }, error = function(e) {
        # 如果出现错误（如 404），处理错误
        if (grepl("404", e$message)) {
            message(paste("File not found (404):", file_name))
        } else {
            message(paste("Error downloading file:", file_name, "Error message:", e$message))
        }
        return(FALSE)  # 返回 FALSE 表示下载失败
    })
}

# 定义提取数据的函数
extract_data_from_vcf <- function(vcf_file, N) {
    # 读取 VCF 文件
    vcf <- read.vcfR(vcf_file)

    # 查看样本列名，确保 sample_name 存在
    sample_columns <- colnames(vcf@gt)
    sample_name <- sample_columns[-1]  # 使用第一个样本列
    print("Available sample columns in the VCF:")
    print(sample_name)  # 打印使用的样本列名

    # 提取固定字段（CHROM, POS, ID, REF, ALT）
    fix_data <- as.data.frame(vcf@fix)

    # 检查 POS 列是否可以正确转换为数值
    print("Checking POS column:")
    print(summary(fix_data$POS))  # 打印 POS 列的摘要信息
    fix_data$POS <- as.numeric(fix_data$POS)
    if (any(is.na(fix_data$POS))) {
        warning("Some positions (POS) could not be converted to numeric. Check the VCF file for inconsistencies.")
    }

    # 提取基因型字段（ES, SE, LP, AF, ID）
    gt_data <- as.data.frame(vcf@gt)

    # 解析 GT 字段中的 ES, SE, LP, AF, ID
    gt_parsed <- do.call(rbind, strsplit(gt_data[, sample_name], ":"))

    # 检查 gt_parsed 的行数是否与 fix_data 的行数一致
    if (nrow(gt_parsed) != nrow(fix_data)) {
        stop("The number of rows in the genotype data does not match the fixed fields. Check the VCF file for inconsistencies.")
    }

    # 创建一个 data.frame，符合 MRlap 的要求
    df <- data.frame(
        chr = fix_data$CHROM,               # 染色体编号
        pos = fix_data$POS,                 # 位置
        rsid = fix_data$ID,                 # SNP 标识符
        ref = fix_data$REF,                 # 参考等位基因
        alt = fix_data$ALT,                 # 替代等位基因
        beta = as.numeric(gt_parsed[, 1]),  # 效应量 (ES)
        se = as.numeric(gt_parsed[, 2]),    # 标准误差 (SE)
        zscore = as.numeric(gt_parsed[, 1]) / as.numeric(gt_parsed[, 2]),  # 计算 Z 分数
        N = N                               # 样本量 (TotalControls + TotalCases)
    )

    print(head(df))  # 查看生成的 data.frame
    print("Data extraction complete.")

    return(df)
}

# 定义基础 URL 和文件名
base_url <- "https://gwas.mrcieu.ac.uk/files/"
exposure_file_name <- "{Exposure_id}.vcf.gz"
outcome_file_name <- "{Outcome_id}.vcf.gz"

# 定义暴露和结果数据的样本量
"""
def MRtool_MRlap(Exposure_id, Outcome_id, path, N_exposure, N_outcome):
    """
    Run MRlap analysis for sample overlap correction using standalone R script.
    Uses subprocess for clean separation and proper error handling.
    """
    import subprocess
    import sys

    # Get the path to the R script relative to this file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path = os.path.join(script_dir, 'r_scripts', 'run_mrlap.R')

    # Build the command
    cmd = [
        'Rscript',
        '--vanilla',
        r_script_path,
        Exposure_id,
        Outcome_id,
        path,
        str(N_exposure),
        str(N_outcome)
    ]

    print(f"[MRtool_MRlap] Running command: {' '.join(cmd)}")

    # Run the process and capture output
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Print output for debugging
    if result.stdout:
        print("[MRtool_MRlap] R stdout:\n", result.stdout)
    if result.stderr:
        print("[MRtool_MRlap] R stderr:\n", result.stderr)

    if result.returncode != 0:
        raise RuntimeError(f"MRlap analysis failed with exit code {result.returncode}. Check R output above.")

    print("[MRtool_MRlap] MRlap analysis completed successfully")


def get_synonyms(term, api_key):
    try:
        # 获取cui
        url = "https://uts-ws.nlm.nih.gov/rest/search/current?apiKey={apiKey}&string={term}&pageNumber=1&pageSize=1".format(
            apiKey=api_key, term=term)
        payload = {}
        headers = {}
        response = requests.request("GET", url, headers=headers, data=payload)
        # print(response.text.encode('utf8'))
        cui = response.json()["result"]["results"][0]["ui"]
        # print(cui)

        # 获取同义词
        url = "https://uts-ws.nlm.nih.gov/rest/content/current/CUI/{cui}/atoms?apiKey={apiKey}&ttys=&language=ENG&pageSize=25".format(
            apiKey=api_key, cui=cui)
        payload = {}
        headers = {}
        response = requests.request("GET", url, headers=headers, data=payload)
        # print(response.json())
        synonyms = [result["name"] for result in response.json()["result"]]

        # 全体小写
        synonyms = [synonym.lower() for synonym in synonyms]

        # 去重
        synonyms = list(set(synonyms))

    except Exception as e:
        print(e)
        synonyms = []

    # print(len(synonyms))
    print(synonyms)

    return synonyms


# ========== GWAS Catalog API functions ==========

def check_keyword_in_gwas_catalog(keyword):
    """Check if a keyword exists in GWAS Catalog via API"""
    url = "https://www.ebi.ac.uk/gwas/rest/api/studies/search"
    params = {
        "pubmedId": "",
        "studyId": "",
        "trait": keyword
    }
    try:
        response = requests.get(url, params=params)
        if response.status_code == 200:
            data = response.json()
            count = data.get("_embedded", {}).get("studies", []).__len__()
            if count > 0:
                print(f"'{keyword}' found in GWAS Catalog ({count} studies)")
                return True
            else:
                print(f"No studies related to '{keyword}' found in GWAS Catalog")
                return False
        else:
            print(f"GWAS Catalog API error: Status code {response.status_code}")
            return None
    except Exception as e:
        print(f"Error querying GWAS Catalog: {e}")
        return None


def get_gwas_id_gwas_catalog(keyword, max_results=30):
    """Get GWAS study information from GWAS Catalog via API"""
    url = "https://www.ebi.ac.uk/gwas/rest/api/studies/search"
    params = {
        "trait": keyword,
        "size": max_results
    }
    results = []
    try:
        response = requests.get(url, params=params)
        if response.status_code == 200:
            data = response.json()
            studies = data.get("_embedded", {}).get("studies", [])
            for study in studies:
                study_info = {
                    "study_id": study.get("studyId"),
                    "trait": study.get("diseaseTrait", {}).get("trait"),
                    "pmid": study.get("pubmedId"),
                    "author": study.get("author"),
                    "publication_date": study.get("publicationDate"),
                    "sample_size": study.get("initialSampleSize"),
                    "platform": study.get("genotypingPlatform")
                }
                results.append(json.dumps(study_info))
            print(f"Found {len(results)} studies in GWAS Catalog for '{keyword}'")
            return results
        else:
            print(f"GWAS Catalog API error: Status code {response.status_code}")
            return []
    except Exception as e:
        print(f"Error fetching from GWAS Catalog: {e}")
        return []


def get_associations_gwas_catalog(study_id):
    """Get SNP associations for a specific study from GWAS Catalog"""
    url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/{study_id}/associations"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()
        else:
            print(f"GWAS Catalog associations API error: {response.status_code}")
            return None
    except Exception as e:
        print(f"Error fetching associations: {e}")
        return None


# ========== FinnGen API functions ==========

def check_keyword_in_fingen(keyword):
    """Check if a keyword exists in FinnGen via API"""
    url = "https://api.finngen.fi/api/phenos"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            phenos = response.json()
            matches = [p for p in phenos if keyword.lower() in p.get("name", "").lower()]
            if len(matches) > 0:
                print(f"'{keyword}' found in FinnGen ({len(matches)} phenotypes)")
                return True
            else:
                print(f"No phenotypes related to '{keyword}' found in FinnGen")
                return False
        else:
            print(f"FinnGen API error: Status code {response.status_code}")
            return None
    except Exception as e:
        print(f"Error querying FinnGen: {e}")
        return None


def get_gwas_id_fingen(keyword, max_results=30):
    """Get phenotype information from FinnGen via API"""
    url = "https://api.finngen.fi/api/phenos"
    results = []
    try:
        response = requests.get(url)
        if response.status_code == 200:
            phenos = response.json()
            matches = [p for p in phenos if keyword.lower() in p.get("name", "").lower()]
            if len(matches) > max_results:
                matches = matches[:max_results]
            for p in matches:
                pheno_info = {
                    "phenocode": p.get("phenocode"),
                    "name": p.get("name"),
                    "description": p.get("description"),
                    "category": p.get("category"),
                    "num_cases": p.get("num_cases"),
                    "population": "Finnish"
                }
                results.append(json.dumps(pheno_info))
            print(f"Found {len(results)} phenotypes in FinnGen for '{keyword}'")
            return results
        else:
            print(f"FinnGen API error: Status code {response.status_code}")
            return []
    except Exception as e:
        print(f"Error fetching from FinnGen: {e}")
        return []


def get_finnen_pheno_info(phenocode):
    """Get detailed information for a specific FinnGen phenotype"""
    url = f"https://api.finngen.fi/api/phenos/{phenocode}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()
        else:
            print(f"FinnGen phenotype API error: {response.status_code}")
            return None
    except Exception as e:
        print(f"Error fetching FinnGen phenotype info: {e}")
        return None


# ========== UK Biobank API functions (via OpenGWAS) ==========

def check_keyword_in_ukbiobank(keyword):
    """Check if a keyword exists in UK Biobank GWAS data via OpenGWAS API"""
    url = "https://gwas-api.mrcieu.ac.uk/search"
    params = {
        "q": keyword,
        "pop": "European"
    }
    try:
        response = requests.get(url, params=params)
        if response.status_code == 200:
            data = response.json()
            # Filter for UK Biobank datasets (usually start with ukb-)
            ukb_results = [d for d in data if d.get("id", "").startswith("ukb-") or "UK Biobank" in d.get("consortium", "")]
            if len(ukb_results) > 0:
                print(f"'{keyword}' found in UK Biobank ({len(ukb_results)} datasets)")
                return True
            else:
                print(f"No UK Biobank datasets found for '{keyword}'")
                return False
        else:
            print(f"UK Biobank (OpenGWAS) API error: Status code {response.status_code}")
            return None
    except Exception as e:
        print(f"Error querying UK Biobank: {e}")
        return None


def get_gwas_id_ukbiobank(keyword, max_results=30):
    """Get UK Biobank GWAS datasets via OpenGWAS API"""
    url = "https://gwas-api.mrcieu.ac.uk/search"
    params = {
        "q": keyword
    }
    results = []
    try:
        response = requests.get(url, params=params)
        if response.status_code == 200:
            data = response.json()
            # Filter for UK Biobank datasets
            ukb_results = [d for d in data if d.get("id", "").startswith("ukb-") or "UK Biobank" in d.get("consortium", "")]
            if len(ukb_results) > max_results:
                ukb_results = ukb_results[:max_results]
            for d in ukb_results:
                dataset_info = {
                    "gwas_id": d.get("id"),
                    "trait": d.get("trait"),
                    "year": d.get("year"),
                    "consortium": d.get("consortium"),
                    "sample_size": d.get("sample_size"),
                    "nsnp": d.get("n snp"),
                    "population": d.get("population")
                }
                results.append(json.dumps(dataset_info))
            print(f"Found {len(results)} UK Biobank datasets for '{keyword}'")
            return results
        else:
            print(f"UK Biobank API error: Status code {response.status_code}")
            return []
    except Exception as e:
        print(f"Error fetching UK Biobank datasets: {e}")
        return []
#可以geneticccccccc
import subprocess
import os

def run_genetic_script(script_name, genetic_dir):
    """
    运行本地遗传学 Shell 脚本，并实时打印日志防卡死
    """
    script_path = os.path.join(genetic_dir, "src", script_name)
    if not os.path.exists(script_path):
        raise FileNotFoundError(f"❌ 找不到遗传学脚本: {script_path}")
        
    cmd = ['bash', script_path]
    print(f"🧬 [Genetic Tool] 执行命令: {' '.join(cmd)}")
    
    try:
        # 使用 Popen 实时捕获长耗时输出
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,cwd=genetic_dir)
        for line in process.stdout:
            print(f"[Genetic Log] {line.strip()}")
        process.wait()
        
        if process.returncode != 0:
            raise RuntimeError(f"脚本 {script_name} 执行失败，退出码 {process.returncode}")
        return True
    except Exception as e:
        print(f"❌ 遗传学流程执行异常: {e}")
        raise e