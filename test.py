import random

parent1 = [
    "ANKS3",
    "FANCE",
    "FANCD2",
    "BLM",
    "KLHL12",
    "IL27",
    "KIAA0430",
    "C16orf3",
    "UFM1",
    "CCDC73",
    "HIATL2",
    "TANC2",
]
parent2 = [
    "NLRC3",
    "NR4A3",
    "ELF5",
    "EXOSC8",
    "MESP2",
    "KIAA0556",
    "CTD-2510F5.6",
    "RNF8",
    "VHL",
    "TCF25",
    "C16orf88",
    "KDM5B",
]

child = []
index = 0
while index < 12:
    print(index)
    whichParent = random.randint(0, 1)

    if whichParent == 0:
        child.append(parent1[index])
    elif whichParent == 1:
        child.append(parent2[index])
    index += 1

print(child)
