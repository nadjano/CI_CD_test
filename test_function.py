import pandas as pd


def concat_strings(a, b):
    return a + b



def test_concat_strings():
    assert concat_strings('hello', 'world') == 'helloworld'
    assert concat_strings('hello', 'world') == 'hello world'



