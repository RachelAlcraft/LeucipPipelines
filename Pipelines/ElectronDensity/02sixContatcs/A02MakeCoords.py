'''
Author: Rachel Alcraft
Date: 9/12/2021
Last update: 10/12/2021
--------------------------------------------------
'''
import pandas as pd

#***************************************************************************************************************
def runMakeCoords(ID,df_geoemtry):
    print('LeucipPipeline:02SixContacts - 1) Creating coords -----------------------------------------------------------')

    results_tup = []
    results_cols = ['SET','ID','pdb_code', 'CMD']
    idxs = df_geoemtry.index
    for idx in idxs:
        pdb_code = df_geoemtry['pdb_code'][idx].lower()
        id1 = df_geoemtry['ID_AB'][idx].lower()
        id2 = df_geoemtry['ID_AC'][idx].lower()
        id3 = df_geoemtry['ID_AD'][idx].lower()
        id4 = df_geoemtry['ID_AE'][idx].lower()
        id5 = df_geoemtry['ID_AF'][idx].lower()

        chain0, aa0, rid0, atom0, coords0 = df_geoemtry['chain'][idx], df_geoemtry['aa'][idx], df_geoemtry['rid'][idx], df_geoemtry['atom'][idx],df_geoemtry['Coords0'][idx]
        chainA, aaA, ridA, atomA, coordsA = df_geoemtry['chainA'][idx], df_geoemtry['aaA'][idx], df_geoemtry['ridA'][idx], df_geoemtry['atomA'][idx],df_geoemtry['CoordsA'][idx]
        chainB, aaB, ridB, atomB,coordsB = df_geoemtry['chainB'][idx], df_geoemtry['aaB'][idx], df_geoemtry['ridB'][idx], df_geoemtry['atomB'][idx],df_geoemtry['CoordsB'][idx]
        chainC, aaC, ridC, atomC, coordsC = df_geoemtry['chainC'][idx], df_geoemtry['aaC'][idx], df_geoemtry['ridC'][idx], df_geoemtry['atomC'][idx], df_geoemtry['CoordsC'][idx]
        chainD, aaD, ridD, atomD, coordsD = df_geoemtry['chainD'][idx], df_geoemtry['aaD'][idx], df_geoemtry['ridD'][idx], df_geoemtry['atomD'][idx], df_geoemtry['CoordsD'][idx]
        chainE, aaE, ridE, atomE, coordsE = df_geoemtry['chainE'][idx], df_geoemtry['aaE'][idx], df_geoemtry['ridE'][idx], df_geoemtry['atomE'][idx], df_geoemtry['CoordsE'][idx]
        chainF, aaF, ridF, atomF, coordsF = df_geoemtry['chainF'][idx], df_geoemtry['aaF'][idx], df_geoemtry['ridF'][idx], df_geoemtry['atomF'][idx], df_geoemtry['CoordsF'][idx]

        if coords0 != '::' and coordsA != '::' and coordsB != '::' and coordsC != '::' and coordsD != '::' and coordsE != '::' and coordsF != '::':
            cx, cy, cz = round(float(coords0.split(':')[0]),5),round(float(coords0.split(':')[1]),5),round(float(coords0.split(':')[2]),5)
            lx,ly,lz = round(float(coordsA.split(':')[0]),5),round(float(coordsA.split(':')[1]),5),round(float(coordsA.split(':')[2]),5)
            px2, py2, pz2 = round(float(coordsB.split(':')[0]),5),round(float(coordsB.split(':')[1]),5),round(float(coordsB.split(':')[2]),5)
            px3, py3, pz3 = round(float(coordsC.split(':')[0]), 5), round(float(coordsC.split(':')[1]), 5), round(float(coordsC.split(':')[2]), 5)
            px4, py4, pz4 = round(float(coordsD.split(':')[0]), 5), round(float(coordsD.split(':')[1]), 5), round(float(coordsD.split(':')[2]), 5)
            px5, py5, pz5 = round(float(coordsE.split(':')[0]), 5), round(float(coordsE.split(':')[1]), 5), round(float(coordsE.split(':')[2]), 5)
            px6, py6, pz6= round(float(coordsF.split(':')[0]), 5), round(float(coordsF.split(':')[1]), 5), round(float(coordsF.split(':')[2]), 5)

            planes = [[id1,px2,py2,pz2],[id2,px3,py3,pz3],[id3,px4,py4,pz4],[id4,px5,py5,pz5],[id5,px6,py6,pz6]]

            for id,px,py,pz in planes:
                #"SLICESFILE|3nir|5|1|0|20_0.1|9.71_-12.376_15.907:0.975_-15.224_7.167:2.52_-13.584_7.06"
                CMD = "SLICESFILE|" + pdb_code  + "|5|1|0|15_0.075|"
                CMD +=  str(cx) + "_" + str(cy) + "_" + str(cz)
                CMD +=  ":" + str(lx) + "_" + str(ly) + "_" + str(lz)
                CMD += ":" + str(px) + "_" + str(py) + "_" + str(pz)
                results_tup.append([ID,id,pdb_code,CMD])

    ### 7) Save the dataframe
    data_for_leuci = pd.DataFrame(results_tup,columns=['SET','ID','pdb_code','CMD'])
    data_for_leuci.to_csv("Csv/" + ID + "_04_LeucipPlus.csv", index=False)
    print("LeucipPipelines:02SixContacts - Saved to", "Csv/" + ID + "_04_LeucipPlus.csv")


#***************************************************************************************************************
