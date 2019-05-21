import matplotlib.pyplot as plt

proteinString = "HHPHHHPH"
coordinates = [[1,1], [2,1], [2,2], [3,2], [4,2], [4,1], [3,1],[3,0]]

# possible = [[0,1],[-1,2],[0,1],[-1,0],[0,1], [1,0]]
# possible = [[0,0], [-1,1], [0,0], [1,1]]
# possible = [[4,1], [3,2]]
possible = []

fig, ax = plt.subplots()
plt.rcParams.update({'font.size': 20})

colorDict = {"P": "b", "H": "r"}
nameDict = {"P": "polair", "H": "hydrofoob"}

ax.grid(False)

if possible:
    for k in range(len(possible)):
        ax.scatter(possible[k][0], possible[k][1], color = "black", alpha = 0.3, s = 100)

    xDotted = [coordinate[0] for coordinate in possible]
    yDotted = [coordinate[1] for coordinate in possible]

    ax.plot(xDotted, yDotted, "k-", zorder = -1, linewidth = 2, alpha = 0.3)

for i in range(len(proteinString)):
    if proteinString[i] == "P":
        ax.scatter(coordinates[i][0], coordinates[i][1], color = "r", s = 100)
    else:
        ax.scatter(coordinates[i][0], coordinates[i][1], color = "b", s = 100)
    plt.text(coordinates[i][0] + 0.08, coordinates[i][1]+ 0.08, (i + 1), fontsize=9)

# ax.scatter(3,2, color = "black", alpha = 0.3, s= 100)
# plt.text(3+0.08,2+0.08,6,fontsize = 9)
#
# ax.scatter(3,3, color = "black", alpha = 0.3, s= 100)
# plt.text(3+0.08,3+0.08,5,fontsize = 9)

xLines = [coordinate[0] for coordinate in coordinates]
yLines = [coordinate[1] for coordinate in coordinates]

ax.plot(xLines, yLines, "-k", zorder=-1)

# ax.scatter(1,1,color=”b”, s = 100, alpha = 0.3)
#
# ax.scatter(-1,1,color=”b”, alpha = 0.3, s = 100)
# ax.scatter(0,2,color=”b”, alpha = 0.3, s = 100)
#
# ax.plot([0,0],[0,1], 'k-', zorder=-1, linewidth = 2)
# ax.plot([0,1],[1,1], 'k--', zorder=-1, linewidth = 2, alpha = 0.3)


ax.axis('equal')

# plt.xlim(-2, 3)
# plt.ylim(-2, 3)

# ax.set_xticks([-2,-1,0,1,2,3])
# ax.set_yticks([-2,-1,0,1,2,3])
ax.axis('off')

plt.show()
