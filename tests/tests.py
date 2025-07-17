#!/usr/bin/env python3

# DIA TMT QUANTIFICATION - TESTS
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


def test1():
    from tmt_chimerys import main
    from tmt_chimerys import __version

    assert int(__version.split(".")[2]) >= 9


def test2():
    from tmt_spectronaut import main
    from tmt_spectronaut import __version

    assert int(__version.split(".")[2]) >= 9
