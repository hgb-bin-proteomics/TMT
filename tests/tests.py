#!/usr/bin/env python3

# DIA TMT QUANTIFICATION - TESTS
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com


def test1():
    from tmt_chimerys import main  # noqa: F401
    from tmt_chimerys import __version

    assert int(__version.split(".")[0]) >= 1


def test2():
    from tmt_chimerys_dda import main  # noqa: F401
    from tmt_chimerys_dda import __version

    assert int(__version.split(".")[0]) >= 1


def test3():
    from tmt_spectronaut import main  # noqa: F401
    from tmt_spectronaut import __version

    assert int(__version.split(".")[0]) >= 1


def test4():
    from tmt_diann import main  # noqa: F401
    from tmt_diann import __version

    assert int(__version.split(".")[0]) >= 1
