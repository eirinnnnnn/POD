import json

def main():
    with open("best.json", "r") as f:
        data = json.load(f)

    perm = data["perm"]
    n = len(perm)

    # Header
    print(f"{n} {n}")

    # Build permutation matrix
    for i in range(n):
        row = ['0'] * n
        row[perm[i]] = '1'
        print("".join(row))

if __name__ == "__main__":
    main()
