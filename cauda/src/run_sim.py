from argparse import ArgumentParser
from sims.utils import calc_kos, interp_fluxes
import numpy as np
import pandas as pd
import psycopg2 as pg

def main():
    parser = ArgumentParser()
    parser.add_argument('--medium_1', '-m1', type=str, required=True)
    parser.add_argument('--medium_2', '-m2', type=str, required=False)
    parser.add_argument('--interp_factor', '-if', type=float, required=False)
    parser.add_argument('--indices', '-i', nargs='+', type=int, required=True)
    parser.add_argument('--rel_abund', '-ra', nargs='+', type=float, required=False)
    args = parser.parse_args()

    # run simulations
    

    # Connect to database
    conn = pg.connect(dbname="cauda", user="postgres")
    cur = conn.cursor()