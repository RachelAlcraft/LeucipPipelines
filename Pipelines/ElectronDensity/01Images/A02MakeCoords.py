'''
Author: Rachel Alcraft
Date: 9/12/2021
Last update: 10/12/2021
--------------------------------------------------
'''
import pandas as pd

#***************************************************************************************************************
def runMakeCoords(ID,df_geoemtry):
    print('LeucipPipeline:01Images - 1) Creating coords -----------------------------------------------------------')

    results_tup = []
    results_cols = ['SET','ID','pdb_code', 'CMD']
    idxs = df_geoemtry.index
    for idx in idxs:
        pdb_code = df_geoemtry['pdb_code'][idx].lower()
        id1 = df_geoemtry['ID'][idx].lower()
        chain1, aa1, rid1, atom1, coords1 = df_geoemtry['chain'][idx], df_geoemtry['aa'][idx], df_geoemtry['rid'][idx], df_geoemtry['atom1'][idx],df_geoemtry['Coords1'][idx]
        chain2, aa2, rid2, atom2, coords2 = df_geoemtry['chain2'][idx], df_geoemtry['aa2'][idx], df_geoemtry['rid2'][idx], df_geoemtry['atom2'][idx],df_geoemtry['Coords2'][idx]
        chain3, aa3, rid3, atom3,coords3 = df_geoemtry['chain3'][idx], df_geoemtry['aa3'][idx], df_geoemtry['rid3'][idx], df_geoemtry['atom3'][idx],df_geoemtry['Coords3'][idx]
        #id1 = str(aa1) + str(rid1) + str(chain1) + "." + str(atom1) + "_" + str(aa2) + str(rid2) + str(chain2) + "." + str(atom2) + "_" + str(aa3) + str(rid3) + str(chain3) + "." + str(atom3)

        if coords1 != '::' and coords2 != '::' and coords3 != '::':
            cx, cy, cz = float(coords1.split(':')[0]),float(coords1.split(':')[1]),float(coords1.split(':')[2])
            lx,ly,lz = float(coords2.split(':')[0]),float(coords2.split(':')[1]),float(coords2.split(':')[2])
            px, py, pz = float(coords3.split(':')[0]),float(coords3.split(':')[1]),float(coords3.split(':')[2])
            cx,cy,cz = round(cx,5),round(cy,5),round(cz,5)
            lx, ly, lz = round(lx, 5), round(ly, 5), round(lz, 5)
            px, py, pz = round(px, 5), round(py, 5), round(pz, 5)
            #"SLICESFILE|3nir|5|1|0|20_0.1|9.71_-12.376_15.907:0.975_-15.224_7.167:2.52_-13.584_7.06"
            CMD = "SLICESFILE|" + pdb_code  + "|5|1|0|15_0.075|"
            CMD +=  str(cx) + "_" + str(cy) + "_" + str(cz)
            CMD +=  ":" + str(lx) + "_" + str(ly) + "_" + str(lz)
            CMD += ":" + str(px) + "_" + str(py) + "_" + str(pz)

            results_tup.append([ID,id1,pdb_code,CMD])

    ### 7) Save the dataframe
    data_for_leuci = pd.DataFrame(results_tup,columns=['SET','ID','pdb_code','CMD'])
    data_for_leuci.to_csv("Csv/" + ID + "_04_LeucipPlus.csv", index=False)
    print("LeucipPipelines:01Images - Saved to", "Csv/" + ID + "_04_LeucipPlus.csv")


#***************************************************************************************************************
