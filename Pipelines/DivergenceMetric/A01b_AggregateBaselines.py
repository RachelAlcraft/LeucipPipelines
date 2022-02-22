import pandas as pd
from LeucipPy import HtmlReportMaker as hrm
from matplotlib import pyplot as plt
import seaborn as sns

dens = [0.5,1.0,2.0,5.0,10.0,50.0]
iters = 500
csv_file = "Csv/A01_Baseline_three"+str(iters)+'_'+str(dens[0])+".csv"
df = pd.read_csv(csv_file)
for i in range(1,len(dens)):
    csv_file = "Csv/A01_Baseline_three" + str(iters) + '_' + str(dens[i]) + ".csv"
    df2 = pd.read_csv(csv_file)
    df = df.append(df2,ignore_index=True)


html_file = "Html/A01b_Baseline_Compare"+str(iters)+".html"
rep = hrm.HtmlReportMaker("Williams Divergence From Trivial: Comparisons",html_file,cols=4)

xys = []
xys.append(['random','stat','samples','density','0.5'])
xys.append(['random','stat','samples','density','1.0'])
xys.append(['random','stat','samples','density','5.0'])
xys.append(['random','stat','samples','density','10.0'])
xys.append(['samples','stat','random','density','0.5'])
xys.append(['samples','stat','random','density','1.0'])
xys.append(['samples','stat','random','density','5.0'])
xys.append(['samples','stat','random','density','10.0'])
xys.append(['random','stat','density','samples','200'])
xys.append(['random','stat','density','samples','1000'])
xys.append(['random','stat','density','samples','2000'])
xys.append(['random','stat','density','samples','10000'])
xys.append(['samples','stat','density','random','0.0'])
xys.append(['samples','stat','density','random','0.2'])
xys.append(['samples','stat','density','random','0.5'])
xys.append(['samples','stat','density','random','10.0'])
for x,y,h,f,v in xys:
    rowdf = df.query(f + '==' + v )
    print(rowdf)
    fig, ax = plt.subplots()
    sns.lineplot(data=rowdf,x=x,y=y,hue=h,palette='tab10')
    plt.title(f + '=' + v)
    plt.legend(title=h,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0,1)
    rep.addPlotOnly(fig, ax)

rep.printReport()