import re

genes = [
    "MAPK14",
    "BDNF",
    "FANCD2OS",
    "NR2F2",
    "PRDM7",
    "NOMO1",
    "PALB2",
    "C16orf96",
    "AAED1",
    "CSNK1A1L",
    "MYBPH",
    "TLK2",
]
# Step 2: Create a regex pattern
pattern = re.compile("|".join(genes))

# Step 3 and 4: Open the file and iterate over each line
count = 0
with open("static/STRING 1.txt", "r") as file:
    for line in file:
        col1, col2, weight = line.strip().split("\t")
        if pattern.fullmatch(col1) and pattern.fullmatch(col2):
            print(line)
            count += 1

print(f"Number of lines where both columns contain a gene from the list: {count}")
