#! /usr/bin/env python3

import pandas as pd 

pm = ["ReversedBasePM", "SemiRigidBasePM", "UGMBasePM", ] 
dimension = list(range(10, 50, 5))
depth = [0.0, 0.1, 0.3, 0.6, 1.0, 5.0, ]
seedSize = list(range(2, 6, 1))

# generate index for bisar results
index = {}
c = 1
for a in pm:
    index[a] = {}
    for b in depth:
        index[a][b] = c
        c += 1
bisar = pd.read_csv("BISAR_S33.csv",  sep=",", lineterminator="\n", )
xColumn = bisar.iloc[:, 0] * 1000.0
pd.set_option("display.max_rows", None, "display.max_columns", None)

def toExcel(content, sheetName, filename, ):
    with pd.ExcelWriter(filename, mode='a', ) as writer:
        content.to_excel(writer, sheet_name=sheetName, index=False, )

# define a function to collect data in each pavement
def collectPMdata(dim, seed, pmName, dep, bisarData=bisar, xSeq=xColumn, bIndex=index):
    output = xColumn.copy() / 1000.0
    for pm_a in pmName:
        dataSet = pd.DataFrame([])
        # combine various depth data
        for dep_a in dep:
            name = "{}{}_{}_{}{}_S33.rpt".format(pm_a, dim, seed, pm_a[0], dep_a, )
            print(name)
            temp0 = pd.read_csv(name, sep="\\s+", )
            xCoord = xColumn.values.tolist()
            for row_i in range(temp0.shape[0]):
                x_i = temp0.iloc[:, 0][row_i]
                judge = (round(x_i, 1) in xCoord)
                if not judge:
                    temp0 = temp0.drop(row_i, )
                else:
                    xCoord.remove(round(x_i, 1))
            temp0 = pd.DataFrame(temp0.to_dict("list"))
            temp1 = bisar.iloc[:, index[pm_a][dep_a]]
            relaE = pd.DataFrame((temp1 - temp0.iloc[:, 1]) / temp1, columns=["RelaE"], )
            dataSet = pd.concat([dataSet, temp0], axis=1, sort=False)
            dataSet = pd.concat([dataSet, temp1], axis=1, sort=False)
            dataSet = pd.concat([dataSet, relaE], axis=1, sort=False)
        # delete redudant column (x coords)
        columnIndex = dataSet.columns[list(range(0, dataSet.shape[1], 4))]
        dataSet = dataSet.drop(columns=columnIndex)
        output = pd.concat([output, dataSet], axis=1, sort=False)
    # print(output)
    return output

# combine by seed
for seed_i in seedSize:
    output = None
    for dim_i in dimension:
        output = pd.concat([output, collectPMdata(dim_i, seed_i, pm, depth, )], \
            axis=1, sort=False)
    toExcel(
        output,
        "Seed{}_overallData".format(seed_i, ), 
        "Seed_overallData_output.xlsx", 
    )

# combine by dimension
for dim_j in dimension:
    output = None
    for seed_j in seedSize:
        output = pd.concat([output, collectPMdata(dim_j, seed_j, pm, depth, )], \
            axis=1, sort=False)
    toExcel(
        output,
        "Dim{}_overallData".format(dim_j, ), 
        "Dim_overallData_output.xlsx", 
    )