from typing import List, NamedTuple, Tuple
import math

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle


class CircleData(NamedTuple):
    c: Tuple[float, float]
    r: float


v = 0 + 0j
# theta = 0
theta = math.pi
u = math.sqrt((v * v.conjugate()).real + 1) * complex(math.cos(theta),
                                                      math.sin(theta))

print(f"(u, v)={u}, {v}")


def get_gen(circ1: CircleData, circ2: CircleData):
    m = np.mat([
        [1, -complex(*circ1.c)],
        [0, 1],
    ])  # type: ignore
    m = np.mat([
        [v.conjugate(), circ1.r * u.conjugate()],
        [u, circ1.r * v],
    ]) * m
    m = np.mat([
        [circ2.r, 0],
        [0, 1],
    ]) * m  # type: ignore
    m = np.mat([
        [1, complex(*circ2.c)],
        [0, 1],
    ]) * m  # type: ignore
    return m


def moebius_on_point(moeb, p: complex):
    return (moeb[0, 0] * p + moeb[0, 1]) / (moeb[1, 0] * p + moeb[1, 1])


def moebius_on_circ(moeb, circ: CircleData):
    P = complex(*circ.c)
    z = P - complex(
        circ.r * circ.r) / (moeb[1, 1] / moeb[1, 0] + P).conjugate()
    Q = moebius_on_point(moeb, z)
    s = abs(Q - moebius_on_point(moeb, P + circ.r))
    return CircleData((Q.real, Q.imag), s)


CIRCS = [
    CircleData((1, 1), 1.0),
    CircleData((-0.6, 0.6), 0.5),
    CircleData((-1, -1), 0.8),
    CircleData((0.9, -0.9), 0.9),
]
gens = [
    get_gen(CIRCS[0], CIRCS[2]),
    get_gen(CIRCS[1], CIRCS[3]),
    get_gen(CIRCS[2], CIRCS[0]),
    get_gen(CIRCS[3], CIRCS[1]),
]
# print(gens)
group = []
inv = [2, 3, 0, 1]
tag: List[int] = []
LEVMAX = 5

fig, ax = plt.subplots()

patches = []
for i, circ in enumerate(CIRCS):
    group.append(gens[i])
    tag.append(i)
    circle = Circle(*circ)
    patches.append(circle)

num = [1, 5]

for lev in range(1, LEVMAX):
    inew = num[lev]
    for iold in range(num[lev - 1] - 1, num[lev]):
        for j in range(4):
            if j == inv[tag[iold]]:
                continue
            group.append(group[iold] * gens[j])
            tag.append(j)
            inew += 1
    num.append(inew)

print(num)
# print(group)

for i in range(num[-1] - 1):
    for j in range(4):
        if j == inv[tag[i]]:
            continue
        newcirc = moebius_on_circ(group[i], CIRCS[j])
        patches.append(Circle(*newcirc))

p = PatchCollection(patches, alpha=0.4)
ax.add_collection(p)
WIDTH = 2
ax.set_xlim(-WIDTH, WIDTH)
ax.set_ylim(-WIDTH, WIDTH)
ax.set_aspect('equal')

plt.show()
