"""Create an English-labelled local analysis file without changing observations.

Usage: python PYTHON/prepare_analysis_data.py --input PATH --output PATH
The input is a restricted analysis-level file, not a substitute for missing
biomarker calibration or geospatial construction code. Never publish its output.
"""
from __future__ import annotations
import argparse
import csv
import datetime as dt
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadstat

CODE_ROOT = Path(__file__).resolve().parents[1]
CJK = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]")
LABELS = {
 "id":"CHARLS participant identifier (restricted)",
 "id_num":"Numeric participant identifier (restricted)",
 "year":"Survey year",
 "city":"Residence city (restricted)",
 "BA_PhenoAge":"Stored study-specific biological age (years; original calibration unavailable)",
 "BAA_PhenoAge":"Stored study-specific biological age acceleration (years; original calibration unavailable)",
 "delta_BAA_PhenoAge":"Change in biological age acceleration from 2011 to 2015 (years)",
 "age":"Chronological age (years)",
 "wave":"CHARLS survey wave",
 "province_code":"Residence province code (restricted)",
 "city_code":"Residence city code (restricted)",
 "sex":"Sex",
 "education_raw":"Original education category",
 "education":"Harmonised education category",
 "hukou_type":"Household registration type",
 "han_ethnicity":"Han ethnicity",
 "marital_status":"Marital status",
 "urban_rural":"Urban or rural residence",
 "hh_consumption_total":"Total household consumption",
 "hh_consumption_percapita":"Household consumption per capita",
 "hh_income_total":"Total household income",
 "retirement":"Retirement status",
 "employed":"Employment status",
 "pension_insurance":"Pension insurance coverage",
 "medical_insurance":"Medical insurance coverage",
 "smoking_current":"Current smoking",
 "alcohol_current":"Current alcohol consumption",
 "metabolic_activity":"Physical activity: total weekly metabolic equivalents",
 "social_participation_score":"Total frequency score across eight social activities",
 "social_participation":"Participation in social activities each month",
 "cognitive_ability":"Cognitive score (0-21; higher is better)",
 "mental_health":"Depressive symptom score (0-30; higher is worse)",
 "social_isolation1":"Social isolation score (0-4)",
 "social_isolation2":"Social isolation score (0-6)",
 "sleep_quality":"Sleep quality: good if reported sleep is 6-8 hours",
 "bmi":"Body mass index (kg/m2)",
 "hospital":"Any hospital admission in the preceding year",
 "outpatient_visits":"Number of outpatient visits in the preceding month",
 "inpatient_visits":"Number of hospital admissions in the preceding year",
 "outpatient_cost_total":"Total outpatient expenditure",
 "adl_difficulty":"Difficulty in any of six activities of daily living",
 "lat":"City latitude (degrees; restricted location information)",
 "TMP":"Temperature (climate adjustment input)",
 "PER":"Precipitation (climate adjustment input)",
 "SNOW":"Snow (climate adjustment input)",
 "WIND":"Wind (climate adjustment input)",
 "DEM":"Elevation (digital elevation model)",
 "bl_tg":"Triglycerides (mg/dL)",
 "bl_hdl":"HDL cholesterol (mg/dL)",
 "bl_ldl":"LDL cholesterol (mg/dL)",
 "bl_hbalc":"Glycated haemoglobin (%)",
 "region_class1":"Macroregion (east, central, west)",
 "region_class2":"North/south region relative to Qinling-Huaihe line",
 "hu_line":"Region relative to Hu Huanyong line",
 "region":"Geographical region",
 "resilience_A":"Standardised reversed frailty index A; higher indicates better health",
 "resilience_B":"Standardised reversed frailty index B; higher indicates better health",
 "BA_PhenoAge_base_2011":"Baseline 2011 stored study-specific biological age (years)",
 "BAA_PhenoAge_base_2011":"Baseline 2011 stored study-specific biological age acceleration (years)",
 "ba_acceleration_A":"Positive biological age acceleration indicator",
 "dataset_type":"Dataset preparation type",
 "imputation_date":"Date recorded by upstream imputation",
 "imputation_method":"Upstream imputation method description",
}
def label(name):
 if name in LABELS: return LABELS[name]
 if name.startswith("traj15_"):
  return "Low/high longitudinal exposure-pattern group: "+name[7:].replace("_"," ")
 if name.startswith("ipw_"):
  return "Inverse probability weight: "+name[4:].replace("_"," ")
 if name.startswith("joint_exp_"):
  return "Four-category joint EVI and ozone exposure group: "+name[10:].replace("_"," ")
 if name.startswith(("EVI","O3","PM25","PM10","NDVI","FVC")):
  bits=name.replace("mean","preceding ").replace("lag","lag ").replace("_"," ")
  units=" (unitless)" if name.startswith(("EVI","NDVI","FVC")) else " (ug/m3)"
  return (bits+units)[:80]
 return name.replace("_"," ").replace("pa ","physical activity ").capitalize()[:80]

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare(src, dst, dictionary_dir):
 src=src.resolve(); dst=dst.resolve()
 if src==dst: raise ValueError("Input and output must be different files.")
 frame, meta=pyreadstat.read_dta(str(src),apply_value_formats=False)
 original=frame.copy(deep=True)
 mappings=json.loads((CODE_ROOT/"DATA"/"source_name_mapping.json").read_text(encoding="ascii"))
 renames={c:mappings[c] for c in frame.columns if c in mappings and c!=mappings[c]}
 if len(set(renames.get(c,c) for c in frame.columns)) != len(frame.columns):
  raise ValueError("Source-name mapping creates duplicate columns; resolve explicitly.")
 frame=frame.rename(columns=renames)
 invalid=[c for c in frame if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]{0,31}",c)]
 if invalid: raise ValueError("Unmapped/non-English source variables: "+", ".join(invalid))
 # Refuse untranslated value labels rather than inventing category meanings.
 values={renames.get(k,k):v for k,v in meta.variable_value_labels.items()}
 unresolved={k:v for k,v in values.items() if any(CJK.search(str(t)) for t in v.values())}
 if unresolved: raise ValueError("Non-English value labels require explicit translation: "+", ".join(unresolved))
 labels={c:label(c) for c in frame}
 dst.parent.mkdir(parents=True,exist_ok=True)
 pyreadstat.write_dta(frame,str(dst),version=15,column_labels=labels,
                     variable_value_labels=values,file_label="Restricted CHARLS analysis input; English metadata")
 check,_=pyreadstat.read_dta(str(dst),apply_value_formats=False)
 pd.testing.assert_frame_equal(frame,check,check_dtype=False,check_exact=True)
 dictionary_dir.mkdir(parents=True,exist_ok=True)
 with (dictionary_dir/"codebook.csv").open("w",encoding="utf-8",newline="") as f:
  writer=csv.writer(f)
  writer.writerow(["variable","label","storage_type","missing_count","value_codes","paper_role"])
  for c in frame:
   role="Supplied outcome" if c in ("BA_PhenoAge","BAA_PhenoAge") else ("Primary exposure" if c in ("EVI_1KM_mean4y","O3_1KM_mean4y","EVI_500m_mean4y","O3_500m_mean4y") else "Identifier/location (restricted)" if c in ("id","id_num","city","lat","city_code","province_code") else "Covariate, sensitivity variable or preparation metadata")
   writer.writerow([c,labels[c],str(frame[c].dtype),int(frame[c].isna().sum()),
                    json.dumps(values.get(c,{}),ensure_ascii=True),role])
 audit={"input_sha256":sha(src),"output_sha256":sha(dst),"rows":len(frame),"columns":len(frame.columns),
        "numeric_values_and_missingness_exact":True,"variable_names_ascii":True,
        "labels_english":not any(CJK.search(x) for x in labels.values()),
        "renamed_variables":renames,"scope":"Metadata/name conversion only; no outcome or estimate recalculation."}
 return audit

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--input",type=Path,required=True)
 p.add_argument("--output",type=Path,required=True)
 p.add_argument("--dictionary-dir",type=Path,default=CODE_ROOT/"DATA"/"analysis_ready")
 p.add_argument("--audit",type=Path)
 a=p.parse_args()
 audit=prepare(a.input,a.output,a.dictionary_dir)
 if a.audit:
  a.audit.parent.mkdir(parents=True,exist_ok=True);a.audit.write_text(json.dumps(audit,indent=2),encoding="utf-8")
 print(json.dumps({k:v for k,v in audit.items() if k!="renamed_variables"},indent=2))
if __name__=="__main__": main()

