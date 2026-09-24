import re

keep = {
    "GSM2719894",
    "GSM2719895",
    "GSM2719896",
    "GSM2719897"
}

gsm=None
title=None
stic=None
srx=None

with open("GSE102094_family.soft") as f:

    for line in f:

        if line.startswith("^SAMPLE"):
            if gsm in keep:
                print(gsm, title, stic, srx, sep="\t")

            gsm=line.split("=")[1].strip()
            title=None
            stic=None
            srx=None

        elif line.startswith("!Sample_title"):
            title=line.split("=",1)[1].strip()

        elif "stic:" in line:
            stic=line.split("stic:")[1].strip()

        elif "Sample_relation" in line and "SRX" in line:
            m=re.search(r"SRX\d+", line)
            if m:
                srx=m.group(0)

    if gsm in keep:
        print(gsm, title, stic, srx, sep="\t")