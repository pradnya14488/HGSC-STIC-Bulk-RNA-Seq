import re

samples = []

with open("GSE102094_family.soft") as f:
    gsm = None
    title = None
    stic = None

    for line in f:

        if line.startswith("^SAMPLE"):
            if gsm:
                samples.append([gsm, title, stic])

            gsm = line.strip().split("=")[1].strip()
            title = None
            stic = None

        elif line.startswith("!Sample_title"):
            title = line.split("=",1)[1].strip()

        elif "stic:" in line:
            stic = line.split("stic:")[1].strip()

    if gsm:
        samples.append([gsm, title, stic])

with open("sample_metadata.tsv","w") as out:
    out.write("GSM\tTITLE\tSTIC\n")

    for s in samples:
        out.write(
            f"{s[0]}\t{s[1]}\t{s[2]}\n"
        )