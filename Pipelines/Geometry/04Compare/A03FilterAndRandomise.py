import pandas as pd


def filterDataframes(set_id,geo,id, df, boundaries):
    print('Filtering dataframes into categories',set_id,id)
    for b in range(0,len(boundaries)):
        ba = boundaries[b]
        if b == 0:
            query='`' + geo +'` <= ' + str(ba)
            print('query=',query)
            df_cut = df.query(query)
            filename = "Csv/" +set_id+ "/" + id + '_' + set_id  + '_' + str(b) + "_03_Geometry.csv"
            print(filename)
            df_cut.to_csv(filename, index=False)
        else:
            bl = boundaries[b-1]
            queryL = '`' + geo + '` > ' + str(bl)
            print('query lower=', queryL)
            queryU = '`' + geo + '` <= ' + str(ba)
            print('query upper =', queryU)
            df_cut = df.query(queryL)
            df_cut = df_cut.query(queryU)
            filename = "Csv/" +set_id+ "/" + id + '_' + set_id  + '_' + str(b) + "_03_Geometry.csv"
            print(filename)
            df_cut.to_csv(filename, index=False)

    bl = boundaries[len(boundaries)-1]
    queryL = '`' + geo + '` > ' + str(bl)
    print('query lower=', queryL)
    df_cut = df.query(queryL)
    filename = "Csv/" +set_id+ "/" + id + '_' + set_id  + '_' + str(len(boundaries)) + "_03_Geometry.csv"
    print(filename)
    df_cut.to_csv(filename, index=False)


def randomise(set_id,id,bound_num):
    concats = []
    print('randomising...',id)
    min_rows = 100000000
    dfs = []
    for i in range(0, bound_num+1):
        filename = "Csv/" + set_id + "/" + id + '_' + set_id + '_' + str(i) + "_03_Geometry.csv"
        df = pd.read_csv(filename)
        rows = len(df.index)
        min_rows = min(min_rows, rows)
        print('Rows=', rows, 'for', filename)
        dfs.append(df)

    print('Create randomised dataframe for', set_id, id, 'with rows=',min_rows)
    for i in range(0, len(dfs)):
        next_df = dfs[i]
        next_df = next_df.sample(n=min_rows)
        concats.append(next_df)

    rand_name = "Csv/" + set_id + '_' + id + "_randomised_04_Geometry.csv"
    df_all = pd.concat(concats,ignore_index=True)
    print('Saving randomised file', rand_name, 'with rows=',len(df_all.index))
    df_all.to_csv(rand_name, index=False)
