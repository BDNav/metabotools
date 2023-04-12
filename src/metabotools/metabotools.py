import pandas as pd

def series2dataframe(s, index_name='index', values_name='values'):
    return pd.DataFrame({index_name: s.index, values_name: s.values})


def get_mrm_name(met, mrm_db):
        mrm_name = mrm_db[mrm_db['MRM Index'] == met]['MRM Name']
        return mrm_name.values.tolist()[0]
        
def db_merge_chem(left, right, how='inner'):
    """Perform a pandas DataFrame merge on 'Chemical Index'.
    
    See pandas.merge for more information.
    
    Parameters
    ----------
    left : DataFrame
    right : DataFrame
    how : {'left', 'right', 'outer', 'inner'}, default 'inner'
        * left: use only keys from left frame (SQL: left outer join)
        * right: use only keys from right frame (SQL: right outer join)
        * outer: use union of keys from both frames (SQL: full outer join)
        * inner: use intersection of keys from both frames (SQL: inner join)
    
    Returns
    -------
    merged : DataFrame
    """
    
    df = pd.merge(left, right, on='Chemical Index', how=how)
    return df



def db_read_excel(excel_file, db_dict):
    """Read excel table(s) into pandas DataFrame(s).
    
    Parameters
    ----------
    excel_file : String of path to excel workbook.
    db_dict : Dictionary with key being database name and value corresponding to the 
        excel sheet name.
        
    Returns
    -------
    result : Dictionary 
        Pandas DataFrame(s) with keys corresponding to input dictionary.
    """
    
    db = {k : pd.read_excel(excel_file, v) for k, v in db_dict.items()}
    return db

def db_merge_chem(left, right, how='inner'):
    """Perform a pandas DataFrame merge on 'Chemical Index'.
    
    See pandas.merge for more information.
    
    Parameters
    ----------
    left : DataFrame
    right : DataFrame
    how : {'left', 'right', 'outer', 'inner'}, default 'inner'
        * left: use only keys from left frame (SQL: left outer join)
        * right: use only keys from right frame (SQL: right outer join)
        * outer: use union of keys from both frames (SQL: full outer join)
        * inner: use intersection of keys from both frames (SQL: inner join)
    
    Returns
    -------
    merged : DataFrame
    """
    
    df = pd.merge(left, right, on='Chemical Index', how=how)
    return df

def db_merge_path(left, right, how='inner'):
    """Perform a pandas DataFrame merge on 'Pathway Index'.
    
    See pandas.merge for more information.
    
    Parameters
    ----------
    left : DataFrame
    right : DataFrame
    how : {'left', 'right', 'outer', 'inner'}, default 'inner'
        * left: use only keys from left frame (SQL: left outer join)
        * right: use only keys from right frame (SQL: right outer join)
        * outer: use union of keys from both frames (SQL: full outer join)
        * inner: use intersection of keys from both frames (SQL: inner join)
    
    Returns
    -------
    merged : DataFrame
    """
    
    df = pd.merge(left, right, on='Pathway Index', how=how)
    return df
    
def db_merge_name(left, right, how='inner'):
    """Perform a pandas DataFrame merge on 'MRM Name'.
    
    See pandas.merge for more information.
    
    Parameters
    ----------
    left : DataFrame
    right : DataFrame
    how : {'left', 'right', 'outer', 'inner'}, default 'inner'
        * left: use only keys from left frame (SQL: left outer join)
        * right: use only keys from right frame (SQL: right outer join)
        * outer: use union of keys from both frames (SQL: full outer join)
        * inner: use intersection of keys from both frames (SQL: inner join)
    
    Returns
    -------
    merged : DataFrame
    """
    
    df = pd.merge(left, right, on='MRM Name', how=how)
    return df

def db_merge_mrm(left, right, how='inner'):
    """Perform a pandas DataFrame merge on 'MRM Index'.
    
    See pandas.merge for more information.
    
    Parameters
    ----------
    left : DataFrame
    right : DataFrame
    how : {'left', 'right', 'outer', 'inner'}, default 'inner'
        * left: use only keys from left frame (SQL: left outer join)
        * right: use only keys from right frame (SQL: right outer join)
        * outer: use union of keys from both frames (SQL: full outer join)
        * inner: use intersection of keys from both frames (SQL: inner join)
    
    Returns
    -------
    merged : DataFrame
    """
    
    df = pd.merge(left, right, on='MRM Index', how=how)
    return df

def db_merge_full(mrm, chemical, pathway):
    """Create a full metabolite database by merging the MRM database, the Chemical 
    database, and the Pathway database.
    
    Parameters
    ----------
    mrm : DataFrame with 'MRM Index' and 'Chemical Index' columns.
    chemical : DataFrame with 'Chemical Index' and 'Pathway Index' columns.
    pathway: DataFrame with 'Pathway Index' column.
    
    Returns
    -------
    merged : DataFrame
    """
    
    df = db_merge_chem(mrm, chemical)
    df = db_merge_path(df, pathway)
    return df
    
def read_excel_sheets(excel_file, sheet_dict):
    return {k : pd.read_excel(excel_file, v) for k, v in sheet_dict.items()}
    
def load_db(db_file):           
    db_dict = {'a': 'Alias',
               'c': 'Chemical',
               'm': 'MRM',
               'p': 'Pathway',
    }

    db_parts = read_excel_sheets(db_file, db_dict)

    df = pd.merge(db_parts['m'], db_parts['c'], on='Chemical Index')
    df = pd.merge(df, db_parts['p'], on='Pathway Index')
    
    return df, db_parts['a']