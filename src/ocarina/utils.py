import numpy as np
from __casac__.table import table as tb

def degreesToRadians(degrees):
    return degrees * np.pi / 180.0

def radiansToDegrees(radians):
    return radians * 180.0 / np.pi

def queryTable(query="", table=""):
    tb_obj = tb()
    tb_obj.open(table)
    query_table = tb_obj.taql(query)
    tb_obj.close()
    return query_table
