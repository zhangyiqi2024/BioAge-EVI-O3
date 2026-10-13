"""Run a portable restricted-data analysis or public aggregate-figure path."""
import argparse,json,os,re,shutil,subprocess,sys,time,platform
from importlib.metadata import version
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUTPUT=ROOT/"output"
LOG_DIR=ROOT/"logs"
def exe(value,names):
 for x in ([value] if value else names):
  p=Path(x)
  if p.exists():return str(p.resolve())
  if shutil.which(x):return shutil.which(x)
 raise RuntimeError("Executable missing; specify --stata / --rscript")
def run(cmd,label,cwd=ROOT):
 env=os.environ.copy()
 if os.name=="nt":env.update({"LANG":"English_United States.utf8","LC_ALL":"English_United States.utf8"})
 env.update({"PYTHONIOENCODING":"utf-8","PYTHONDONTWRITEBYTECODE":"1"})
 startup=None
 if os.name=="nt":
  startup=subprocess.STARTUPINFO();startup.dwFlags|=subprocess.STARTF_USESHOWWINDOW;startup.wShowWindow=0
 print("[RUN]",label,flush=True)
 with (LOG_DIR/(label+".log")).open("w",encoding="utf-8") as f:
  r=subprocess.run(cmd,cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT,startupinfo=startup)
 if r.returncode:raise RuntimeError(label+" failed; see logs/"+label+".log")
 return {"step":label,"returncode":0}
def main():
 global OUTPUT,LOG_DIR
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--mode",choices=("analysis","public","check"),default="public")
 p.add_argument("--output-dir",type=Path,help="Public-mode outputs and logs; defaults to 04_code/output.")
 p.add_argument("--stata");p.add_argument("--rscript")
 p.add_argument("--skip-stata",action="store_true",help="Resume R only after a verified Stata run.")
 a=p.parse_args();start=time.monotonic();history=[]
 if a.output_dir and a.mode!="public":raise ValueError("--output-dir is only supported in public mode")
 if a.mode=="public":
  OUTPUT=a.output_dir.resolve() if a.output_dir else ROOT/"output"
  LOG_DIR=OUTPUT/"logs"
 LOG_DIR.mkdir(parents=True,exist_ok=True)
 if a.mode=="public":
  session={"python":platform.python_version(),"platform":platform.platform(),
           "packages":{name:version(name) for name in ["matplotlib","pandas","numpy","Pillow"]}}
  (LOG_DIR/"session_info.json").write_text(json.dumps(session,indent=2),encoding="utf-8")
  (LOG_DIR/"pipeline_result.json").write_text(json.dumps({"status":"RUNNING","mode":a.mode}),encoding="utf-8")
 if a.mode=="public":
  main_out=OUTPUT/"main"
  si_out=OUTPUT/"supplementary"
  table_out=OUTPUT/"tables"
  history.append(run([sys.executable,str(ROOT/"PYTHON"/"compose_main_figures.py"),"--root",str(ROOT.parent),"--out",str(main_out)],"public_figures",ROOT.parent))
  history.append(run([sys.executable,str(ROOT/"PYTHON"/"render_supplementary_figures.py"),"--root",str(ROOT.parent),"--output-dir",str(si_out)],"public_supplementary_figures",ROOT.parent))
  history.append(run([sys.executable,str(ROOT/"PYTHON"/"render_public_additional.py"),"--root",str(ROOT.parent),"--main-out",str(main_out),"--si-out",str(si_out)],"public_additional_figures",ROOT.parent))
  history.append(run([sys.executable,str(ROOT/"PYTHON"/"render_public_s5.py"),"--root",str(ROOT.parent),"--output-dir",str(si_out)],"public_s5",ROOT.parent))
  history.append(run([sys.executable,str(ROOT/"PYTHON"/"render_public_tables.py"),"--root",str(ROOT.parent),"--output-dir",str(table_out)],"public_tables",ROOT.parent))
  history.append(run([sys.executable,str(ROOT/"PYTHON"/"audit_public_reproduction.py"),"--root",str(ROOT.parent),"--main-output",str(main_out),"--si-output",str(si_out),"--table-output",str(table_out),"--output",str(LOG_DIR/"public_reproduction_status.json")],"public_reproduction_audit",ROOT.parent))
 else:
  r=exe(a.rscript,["Rscript","Rscript.exe"])
  history.append(run([r,"--vanilla","CODE/00_check_environment.R"],"dependency_preflight"))
  if a.mode=="analysis":
   if not (ROOT/"DATA/analysis_ready/analysis_dataset.dta").exists():raise RuntimeError("Restricted analysis input absent; see DATA/BUILD_SCOPE.md.")
   if not a.skip_stata:
    stata=exe(a.stata,["StataMP-64","StataSE-64","stata-mp","stata"])
    wrapper=ROOT/"logs"/"run_stata_wrapper.do"
    wrapper.write_text("""clear all
set more off
capture log close _all
log using "logs/stata_pipeline.log", name(pipeline_wrapper) text replace
capture noisily do "CODE/00_master_run.do"
local result = _rc
display "PIPELINE_RETURN_CODE=%result'"
file open status using "logs/stata_exit_status.txt", write replace text
file write status "%result'" _n
file close status
capture log close _all
exit %result', clear
""".replace("%",chr(96)),encoding="utf-8")
    history.append(run([stata,"/e" if os.name=="nt" else "-b","do",str(wrapper)],"stata_launch"))
    marker=ROOT/"logs"/"stata_exit_status.txt"
    if not marker.exists() or marker.read_text().strip()!="0":raise RuntimeError("Stata failed; inspect logs/stata_pipeline.log and per-step logs.")
   history.append(run([r,"--vanilla","CODE/00_master_run.R"],"r_pipeline"))
 result={"mode":a.mode,"status":"PASS","steps":history,"elapsed_seconds":round(time.monotonic()-start,2)}
 (LOG_DIR/"pipeline_result.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
 print(json.dumps(result,indent=2))
if __name__=="__main__":
 try:main()
 except Exception as e:
  LOG_DIR.mkdir(parents=True,exist_ok=True)
  (LOG_DIR/"pipeline_result.json").write_text(json.dumps({"status":"FAIL","error":str(e)},indent=2),encoding="utf-8")
  print("ERROR:",e,file=sys.stderr);sys.exit(1)

