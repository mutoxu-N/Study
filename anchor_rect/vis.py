import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import *
from PIL import Image
import os

FILE = "sample.txt"
N = 100

SHOW = True
TEXT = False

colors = [
    "#faf",
    "#caf",
    "#aaf",
    "#acf",
    "#aff",
    "#afc",
    "#afa",
    "#cfa",
    "#ffa",
    "#fca",
    "#faa",
    "#fac",
]


ps = None
rs = None
order = list(range(N))

with open(FILE, encoding="utf-8") as f:
    ps = list(
        map(lambda l: list(map(float, l.split())), f.readlines()))

with open(f"./results/{FILE}", "r") as f:
    lines = f.readlines()
    rs = list(
        map(lambda l: list(map(float, l.split())), lines[:N]))

fig = plt.figure(figsize=(7.8, 8))
ax = fig.add_subplot(1, 1, 1)

DELTA = 0.01
ax.set_xlim(-DELTA, 1+DELTA)
ax.set_ylim(-DELTA, 1+DELTA)
ax.set_aspect("equal")

area = .0
im_cnt = 0
ax.add_patch(Rectangle(
    xy=(0, 0),
    width=1,
    height=1,
    color="#FF0000",
    fill=None,
))

remain = set(range(N))
for i in range(len(order)):
    p = ps[order[i]]
    print(p)
    ax.plot(
        p[0],
        p[1],
        ".",
        c=colors[i % len(colors)],
        zorder=1,
    )
    remain.remove(order[i])

for r in remain:
    p = ps[r]
    i += 1
    ax.plot(
        p[0],
        p[1],
        ".",
        c=colors[i % len(colors)],
        zorder=1,
    )

if TEXT:
    for i, p in enumerate(ps):
        ax.text(
            p[0],
            p[1],
            str(i),
            c="#000000",
            zorder=2,
        )


for i in range(len(order)):
    r = rs[order[i]]
    area += (r[2]-r[0])*(r[3]-r[1])
    ax.add_patch(Rectangle(
        xy=r[:2],
        width=(r[2]-r[0]),
        height=(r[3]-r[1]),
        color=colors[i % len(colors)],
        fill=True,
        zorder=0,
    ))
    ax.add_patch(Rectangle(
        xy=r[:2],
        width=(r[2]-r[0]),
        height=(r[3]-r[1]),
        color="#000000",
        fill=None,
        zorder=0,
    ))


ax.set_title(f"{FILE}, area={area:.8}")
plt.savefig(f"./results/{FILE}.png")

if SHOW:
    plt.show()
plt.close()

print("END")
