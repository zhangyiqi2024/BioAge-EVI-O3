"""Export a reversible paper-facing schema; microdata remain local and restricted."""
import argparse,csv,json,re
from pathlib import Path
import pandas as pd
import pyreadstat
ROOT=Path(__file__).resolve().parents[1]
ALIASES={
 "id":"participant_id","id_num":"participant_id_numeric","year":"survey_year","wave":"survey_wave",
 "age":"chronological_age","BAA_PhenoAge":"stored_baa","BA_PhenoAge":"stored_biological_age",
 "BAA_PhenoAge_base_2011":"stored_baa_baseline_2011",
 "BA_PhenoAge_base_2011":"stored_bioage_baseline_2011",
 "delta_BAA_PhenoAge":"change_stored_baa_2011_2015",
 "EVI_1KM_mean4y":"residential_evi_1km_mean4y",
 "EVI_500m_mean4y":"residential_evi_500m_mean4y",
 "O3_1KM_mean4y":"ozone_1km_mean4y",
 "O3_500m_mean4y":"ozone_500m_mean4y"}
def canonical(name):
 return ALIASES.get(name,name.lower())
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--input",type=Path,default=ROOT/"DATA/analysis_ready/analysis_dataset.dta")
 p.add_argument("--output-dir",type=Path,default=ROOT/"DATA/canonical")
 p.add_argument("--restore-internal",type=Path,help="Restore internal names from this canonical DTA.")
 p.add_argument("--restored-output",type=Path)
 a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
 if a.restore_internal:
  if not a.restored_output:raise ValueError("--restored-output is required")
  if a.restored_output.resolve()==a.restore_internal.resolve():raise ValueError("Source and destination must differ")
  frame,meta=pyreadstat.read_dta(str(a.restore_internal))
  with (ROOT/"DATA/canonical_schema_crosswalk.csv").open(encoding="utf-8-sig") as f:rows=list(csv.DictReader(f))
  inverse={r["canonical_variable"]:r["internal_variable"] for r in rows}
  frame=frame.rename(columns=inverse)
  labels={inverse.get(k,k):v for k,v in meta.column_names_to_labels.items()}
  vals={inverse.get(k,k):v for k,v in meta.variable_value_labels.items()}
  pyreadstat.write_dta(frame,str(a.restored_output),version=15,column_labels=labels,variable_value_labels=vals)
  return
 frame,meta=pyreadstat.read_dta(str(a.input),apply_value_formats=False)
 mapping={c:canonical(c) for c in frame}
 if len(set(mapping.values()))!=len(mapping):raise ValueError("Canonical aliases are not one-to-one")
 if any(not re.fullmatch(r"[a-z][a-z0-9_]{0,31}",v) for v in mapping.values()):
  raise ValueError("Canonical variable violates Stata name format")
 canon=frame.rename(columns=mapping)
 labels={mapping[k]:v for k,v in meta.column_names_to_labels.items()}
 labels["stored_baa"]="Stored biological age acceleration (years; supplied study-specific measure)"
 labels["stored_biological_age"]="Stored biological age (years; supplied study-specific measure)"
 vals={mapping[k]:v for k,v in meta.variable_value_labels.items()}
 dest=a.output_dir/"analysis_canonical.dta"
 pyreadstat.write_dta(canon,str(dest),version=15,column_labels=labels,variable_value_labels=vals,
  file_label="Restricted CHARLS microdata; paper-facing canonical English schema")
 # 17 significant digits preserve binary double values on round-trip.
 canon.to_csv(a.output_dir/"analysis_canonical.csv",index=False,float_format="%.17g")
 back,_=pyreadstat.read_dta(str(dest),apply_value_formats=False)
 back=back.rename(columns={v:k for k,v in mapping.items()})
 pd.testing.assert_frame_equal(frame,back,check_dtype=False,check_exact=True)
 with (ROOT/"DATA/canonical_schema_crosswalk.csv").open("w",encoding="utf-8",newline="") as f:
  w=csv.writer(f);w.writerow(["internal_variable","canonical_variable","paper_term","buffer","time_window","label"])
  for c in frame:
   buffer="1 km" if "1km" in mapping[c] else "500 m" if "500m" in mapping[c] else ""
   window="Mean of four preceding years" if "mean4y" in c else ""
   term={"BAA_PhenoAge":"Stored BAA","BA_PhenoAge":"Stored biological age","age":"Chronological age"}.get(c,
         "Residential greenness (EVI)" if c.startswith("EVI") else "Ozone (O3)" if c.startswith("O3") else labels[mapping[c]])
   w.writerow([c,mapping[c],term,buffer,window,labels[mapping[c]]])
 result={"rows":len(frame),"columns":len(frame.columns),"mapping_one_to_one":True,
         "restored_values_and_missingness_exact":True,"public_distribution_allowed":False}
 (a.output_dir/"conversion_verification.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
 print(json.dumps(result))
if __name__=="__main__":main()

