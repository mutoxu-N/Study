import random

N = 100

ps = {(0, 0)}
cnt = 1
while cnt < N:
    x = random.random()
    y = random.random()
    if (x, y) not in ps:
        ps.add((x, y))
        cnt += 1

with open("sample.txt", "w") as f:
    for i, p in enumerate( ps ):
        f.write(f"{p[0]} {p[1]}\n")
