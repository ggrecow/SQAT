"""Joins the data.json of several runs of sqat_report_build (one per folder of
test/, as the parallel jobs of the CI run them) into one report page.

    python3 merge_reports.py out/index.html gui/data.json reference/data.json unit/data.json

A folder whose results were reused from an earlier run (the reference block,
when the code it tests did not change) carries a file named cached next to its
data.json; the page then says so under Notes.
"""
import json
import os
import sys


def block_of(path):
    # the folder of a data.json, without the report- of the CI artifacts
    name = os.path.basename(os.path.dirname(os.path.abspath(path)))
    return name[len("report-"):] if name.startswith("report-") else name


def main(target, sources):
    here = os.path.dirname(os.path.abspath(__file__))
    runs = []
    notes = []
    for path in sources:
        if not os.path.isfile(path):                 # a job that ended before its report
            notes.append(f"{block_of(path)}: no results (the job ended early, see its log)")
            continue
        with open(path, encoding="utf-8") as f:
            runs.append((path, json.load(f)))
    meta = dict(runs[0][1]["meta"])
    for path, d in runs:
        if os.path.isfile(os.path.join(os.path.dirname(path), "cached")):
            m = d["meta"]
            notes.append(f"{block_of(path)}: reused from the run of {m['date']} ({m['git']}), "
                         "since the code it tests did not change")
    if notes:
        meta["notes"] = "; ".join(notes)
    tests = [t for _, d in runs for t in (d["tests"] if isinstance(d["tests"], list) else [d["tests"]])]
    records = "".join(d["records"] for _, d in runs)
    log = "\n".join(f"===== {block_of(p)} =====\n{d['log']}" for p, d in runs)
    data = {"meta": meta, "tests": tests, "records": records, "log": log}
    text = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    with open(os.path.join(here, "template.html"), encoding="utf-8") as f:
        page = f.read().replace("/*__DATA__*/", f"window.REPORT = {text};")
    os.makedirs(os.path.dirname(os.path.abspath(target)), exist_ok=True)
    with open(target, "w", encoding="utf-8") as f:
        f.write(page)
    failed = sum(t["status"] == "failed" for t in tests)
    ignored = sum(t["status"] == "ignored" for t in tests)
    print(f"report: {target} | tests: {len(tests)}, failed {failed}, ignored {ignored}")
    return failed == 0 and len(runs) == len(sources)


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1], sys.argv[2:]) else 1)
