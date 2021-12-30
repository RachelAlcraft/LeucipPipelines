from datetime import datetime as dt

def getTime():
    now = dt.now()
    tm =  now.strftime("%a %d%h%y-%H:%M::")
    return tm