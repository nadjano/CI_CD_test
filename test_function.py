#!/usr/bin/python3
import pandas as pd

"""Module providing a function printing python version."""

def concat_strings(a, b):
    """
    Concatenate two strings
    """
    return a + b

def test_concat_strings():
    """
    Test the concat_strings function
    """
    assert concat_strings('hello', 'world') == 'helloworld'
    assert concat_strings('hello', 'world') == 'hello world'
