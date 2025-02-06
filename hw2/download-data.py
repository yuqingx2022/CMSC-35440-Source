import json
import os
import re
import shutil
import tarfile
from pathlib import Path

import pandas as pd
import requests
from tqdm import trange

max_size = 2000
batch_size = 10
files_endpt = "https://api.gdc.cancer.gov/files"
data_endpt = "https://api.gdc.cancer.gov/data"
suffix = "rna_seq.augmented_star_gene_counts.tsv"

os.makedirs("data", exist_ok=True)

filters = {
    "op": "and",
    "content": [
        {
            "op": "in",
            "content": {
                "field": "cases.project.project_id",
                "value": ["TCGA-LUAD", "TCGA-LUSC"],
            },
        },
        {
            "op": "in",
            "content": {
                "field": "files.experimental_strategy",
                "value": ["RNA-Seq"],
            },
        },
        {
            "op": "in",
            "content": {
                "field": "files.data_format",
                "value": ["tsv"],
            },
        },
        {
            "op": "in",
            "content": {
                "field": "files.data_category",
                "value": ["transcriptome profiling"],
            },
        },
        {
            "op": "in",
            "content": {
                "field": "files.data_type",
                "value": ["Gene Expression Quantification"],
            },
        },
        {
            "op": "in",
            "content": {
                "field": "files.access",
                "value": ["open"],
            },
        },
        {
            "op": "in",
            "content": {
                "field": "cases.samples.tissue_type",
                "value": ["tumor"],
            },
        },
        {
            "op": "in",
            "content": {
                "field": "cases.samples.sample_type",
                "value": ["primary tumor"],
            },
        },
    ],
}

params = {
    "filters": json.dumps(filters),
    "fields": ",".join(
        [
            "file_id",
            "file_name",
            "cases.submitter_id",
            "cases.samples.submitter_id",
            "cases.samples.portions.submitter_id",
            "cases.samples.portions.analytes.submitter_id",
            "cases.samples.portions.analytes.aliquots.submitter_id",
            "cases.project.project_id",
        ]
    ),
    "format": "JSON",
    "size": max_size,
}

response = requests.get(files_endpt, params=params)

data = []
for x in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
    case = x["cases"][0]
    sample = case["samples"][0]
    portion = sample["portions"][0]
    analyte = portion["analytes"][0]
    aliquot = analyte["aliquots"][0]

    data.append(
        {
            "file_id": x["file_id"],
            "file_name": x["file_name"],
            "case_id": case["submitter_id"],
            "sample_id": sample["submitter_id"],
            "portion_id": portion["submitter_id"],
            "analyte_id": analyte["submitter_id"],
            "aliquots_id": aliquot["submitter_id"],
            "project_id": case["project"]["project_id"],
        }
    )
df = pd.DataFrame(data)
assert len(df) < max_size
assert df["file_name"].str.endswith(suffix).all()

df.to_csv("data/metadata.csv", index=False)

print(
    f"Average number of sequencing files per case:",
    f"{df.groupby('case_id').size().mean():0.3f}",
)

for i in trange(0, len(df), batch_size):
    file_ids = df["file_id"].iloc[i : i + batch_size].to_list()
    params = {"ids": file_ids}

    response = requests.post(
        data_endpt,
        data=json.dumps(params),
        headers={"Content-Type": "application/json"},
    )

    if len(file_ids) == 1:
        response_head_cd = response.headers["Content-Disposition"]
        file_name = re.findall("filename=(.+)", response_head_cd)[0]
        if file_name.endswith(suffix):
            with open("data/" + file_name, "wb") as output_file:
                output_file.write(response.content)
    else:
        with open("data/temp.tar.gz", "wb") as output_file:
            output_file.write(response.content)
        with tarfile.open("data/temp.tar.gz", "r") as tar:
            tar.extractall("data/temp")

        for root, dirs, files in os.walk("data/temp"):
            root = Path(root)
            for f in files:
                if f.endswith(suffix):
                    shutil.move(root / f, Path("data") / f)

        shutil.rmtree("data/temp")
        os.remove("data/temp.tar.gz")
