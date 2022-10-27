import os
from casatools import table

tb = table()


def query_table(query: str = None, table_name: str = None):
    if query is None and table_name is None:
        raise ValueError("Query or table cannot be None")
    else:
        if os.path.exists(table_name):
            tb.open(table_name)
        else:
            raise ValueError("Table name does not exist")
        result = tb.taql(query)
        tb.close()
        return result
