from __future__ import print_function, division

import sys
import time
import subprocess
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from astropy.time import Time

sys.path.append('/home/rgm/soft/python/lib/')
from librgm.xmatch import match_lists
from librgm.plotid import plotid

def mk_corner_file(file, band, exptime = False):
    """

    """
    t = Table.read(file)
    lines = "#Corner file for VHS from " + file + "\n#RA_min    RA_max    Dec_min    Dec_max    Framesetid\n"
    n = 0

    if band + "MINRA" in t.columns:
        while n < len(t):
            if t[band + "MINRA"][n] > -999999000.0:
                ra_min = str(t[band + "MINRA"][n])
                ra_max = str(t[band + "MAXRA"][n])
                dec_min = str(t[band + "MINDEC"][n])
                dec_max = str(t[band + "MAXDEC"][n])
                id = str(t["FRAMESETID"][n])
                lines += (ra_min + "    " + ra_max + "    " + dec_min + "    " + dec_max + "    " + id + "\n")
            n += 1

    else:
        if band == "KS":
            band = "Ks"
        if not exptime:
            ids = np.where( (t["filtname"] == band) )[0]
        else:
            print(exptime)
            ids = np.where( (t["filtname"] == band) & (t["image_exptime"] == exptime) )[0]
        if band == "all":
            ids = np.ones(len(t), dtype = np.bool)
        t = t[ids]
        while n < len(t):
            ra_max = str(max([t["ra1"][n], t["ra2"][n], t["ra3"][n], t["ra4"][n]]))
            ra_min = str(min([t["ra1"][n], t["ra2"][n], t["ra3"][n], t["ra4"][n]]))
            if float(ra_max) - float(ra_min) > 10.0:
                ra_min, ra_max = ra_max, ra_min
            dec_max = str(max([t["dec1"][n], t["dec2"][n], t["dec3"][n], t["dec4"][n]]))
            dec_min = str(min([t["dec1"][n], t["dec2"][n], t["dec3"][n], t["dec4"][n]]))
            id = str(t["image_id"][n])
            lines += (ra_min + "    " + ra_max + "    " + dec_min + "    " + dec_max + "    " + id + "\n")

            n += 1
    try:
        with open(file[:-5] + "_corners_" + band + ".txt", "w") as f:
            f.write(lines)
        print("Made corner file:", file[:-5] + "_corners_" + band + ".txt")
        return file[:-5] + "_corners_" + band + ".txt"
    except IOError:
        with open(file[:9] + file[18:-5] + "_corners_" + band + ".txt", "w") as f:
            f.write(lines)
        print("Made corner file:", file[:9] + file[18:-5] + "_corners_" + band + ".txt")
        return file[:9] + file[18:-5] + "_corners_" + band + ".txt"

def mk_corner_file_all(file):

    lines = ""

    c = 0
    print(file)
    for band in ["Y", "J", "H", "KS"]:
        print(file)
        print(file[:-5])
        # 20180504_rgm: added cfile
        cfile = file[:-5] + "_corners_" + band + ".txt"
        print('open:', cfile)
        f = open(cfile, "r")

        if c == 0:
            for line in f:
                lines += (line + "\n")
            c += 1
        else:
            for line in f:
                if line[0] <> "#":
                    lines += (line + "\n")

    with open(file[:-5] + "_corners_all.txt", "w") as f:
        f.write(lines)

    return file[:-5] + "_corners_all.txt"

def mk_weight_file(file, band, conv):

    t = Table.read(file)
    lines = "#Weight file for VHS from " + file + "\n#Weight    Framesetid\n"
    n = 0
    if band + "MINRA" in t.columns:
            while n < len(t):
                if t[band + "MINRA"][n] > -999999000.0:
                    weight = str(t[band + "_DEPTH_DYE2006"][n])
                    id = str(t["FRAMESETID"][n])
                    lines += (weight + "    " + id + "\n")
                n += 1
    else:
        if band == "KS":
            band = "Ks"
        if not exptime:
            ids = np.where( (t["filtname"] == band) )[0]
        else:
            print(exptime)
            ids = np.where( (t["filtname"] == band) & (t["image_exptime"] == exptime) )[0]
        if band == "all":
            ids = np.ones(len(t), dtype = np.bool)
        t = t[ids]
        while n < len(t):
            weight = str(t["maglim"][n])
            id = str(t["image_id"][n])
            lines += (weight + "    " + id + "\n")
            n += 1

    with open(file[:-5] + "_weights_" + band + ".txt", "w") as f:
        f.write(lines)

    print("Made weight file:", file[:-5] + "_weights_" + band + ".txt")
    return file[:-5] + "_weights_" + band + ".txt"

def mk_weight_file_all(file):

    lines = ""

    c = 0
    print(file)
    for band in ["Y", "J", "H", "KS"]:
        f = open(file[:-5] + "_weights_" + band + ".txt", "r")

        if c == 0:
            for line in f:
                lines += (line + "\n")
            c += 1
        else:
            for line in f:
                if line[0] <> "#":
                    lines += (line + "\n")

    with open(file[:-5] + "_weights_all.txt", "w") as f:
        f.write(lines)

    return file[:-5] + "_weights_all.txt"

def run_mangle(file, band, conv, exptime = False):
    """


    """

    print('Running mangle on:', file)
    if band == "all":
        corner_file = mk_corner_file_all(file)
    else:
        corner_file = mk_corner_file(file, band, exptime = exptime)

    corners_out = corner_file[:-4] + "_mangle_output_" + band
    print(corner_file)
    print(corners_out)
    subprocess.call([manglepath + "poly2poly",
                     "-ir", corner_file, corners_out])

    if band == "all":
        weight_file = mk_weight_file_all(file)
    else:
        weight_file = mk_weight_file(file, band, conv)
    weights_out = weight_file[:-4] + "_mangle_output_" + band
    subprocess.call([manglepath + "weight",
                     "-z", weight_file, corners_out, weights_out])

    pix_out = file[:-5] + "_pix_mangle_out_" + band
    subprocess.call([manglepath + "pixelize",
                     weights_out, pix_out])

    snap_out = file[:-5] + "_snap_mangle_out_" + band
    subprocess.call([manglepath + "snap",
                     pix_out, snap_out])

    balk_out = file[:-5] + "_balk_mangle_out_" + band
    subprocess.call([manglepath + "balkanize",
                     snap_out, balk_out, "-Bx"])

    unify_out = file[:-5] + "_unify_mangle_out_" + band
    subprocess.call([manglepath + "unify",
                     balk_out, unify_out])

    area_file = file[:-5] + "_area_mangle_out_" + band
    subprocess.call([manglepath + "poly2poly",
                     "-oa", unify_out, area_file])

    depth_file = file[:-5] + "_depth_mangle_out_" + band
    subprocess.call([manglepath + "poly2poly",
                     "-ow", unify_out, depth_file])

    return area_file, depth_file

def area_depth(file, area_file, depth_file, band, conv):

    areas = []
    print(area_file)
    with open(area_file, "r") as af:
        for line in af:
            area = line.split()[0]
            if area[0] == "0":
                areas.append(float(area))

    depths = []
    with open(depth_file, "r") as df:
        for line in df:
            depth = line.split()[0]
            if depth[0] <> "w":
                depths.append(float(depth))

    areas = np.array(areas)*((180.0/np.pi)**2)
    depths = np.array(depths) + conv

    print("Total area:", np.sum(areas))

    bins = np.linspace(min(depths), max(depths), 101)
    areas_out = np.zeros(100)
    n = 0
    while n < len(bins)-1:
        ids = np.where( (depths >= bins[n]) & (depths < bins[n + 1]) )[0]
        areas_out[n] += np.sum(areas[ids])
        n += 1

    cum_areas = np.zeros(100)
    n = 0
    for area in areas_out[::-1]:
        cum_areas[n] = area + cum_areas[n-1]
        n += 1

    tot_area = np.sum(areas)
    half_area = tot_area/2.0
    n = 0
    while n < len(cum_areas)-1:
        if half_area > cum_areas[n] and half_area < cum_areas[n+1]:
            half_lim = bins[::-1][1:][n]
        n += 1
    print("Magnitude limit for half the area:", half_lim)

    half_lim = "%0.2f" % half_lim
    tot_area = "%0.2f" % tot_area

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(bins[::-1][1:], cum_areas)
    ax.set_title("Area in the " + band + " band\n from " + file)
    print(conv)
    if conv <= 0.0:
        ax.set_xlabel("Magnitude Limit " + band + " band [Vega]")
    elif conv > 0.0:
        ax.set_xlabel("Magnitude Limit " + band + " band [AB]")
    ax.set_ylabel("Cumulative Area / sq degrees")
    ax.text(0.1, 0.6, "Total Area: " + tot_area + r" deg$^2$" + "\nMag lim for 50%: " + half_lim, transform = ax.transAxes, bbox = dict(facecolor = "white", alpha = 0.7))
    plt.savefig(file[:-5] + "_area_depth_" + band + ".png")
    plt.show()

def area_depth_all(file, system = "AB", norm = False, save_data = False):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    bands = ["Y", "J", "H", "Ks"]
    #bands = ["J", "H", "Ks"]
    mag_lim_list = []
    colours = ["b", "g", "y", "r"]
    #colours = ["orange", "red", "darkred"]
    info = "     Area  Mag Lim 50%\n"
    m = 0
    data = []
    #if file == "/data/vhs/vsa/dqc/VHSv20140517/vhs_vsa_dqc_tiles_fs_metadata.fits":
    if "vsa" in file:
        if system == "AB":
            convs = [0.618-0.55, 0.937-0.55, 1.384-0.55, 1.839-0.55]
        else:
            convs = [-0.55, -0.55, -0.55, -0.55]
    else:
        if system == "AB":
            convs = [0.618, 0.937, 1.384, 1.839]
        else:
            convs = [0.0, 0.0, 0.0, 0.0]

    t = Table.read(file)

    if "YMINRA" in t.columns:
        min_mjd = min(min(t["YMJDOBS"]), min(t["JMJDOBS"]), min(t["HMJDOBS"]), min(t["KSMJDOBS"]))
        max_mjd = max(max(t["YMJDOBS"]), max(t["JMJDOBS"]), max(t["HMJDOBS"]), max(t["KSMJDOBS"]))
        t_mjd = Time([min_mjd, max_mjd], scale = "utc", format = "mjd")
        mjd_info = "MJD Range: %0.2f:%0.2f\n" % (min_mjd, max_mjd)
        mjd_info += (t_mjd.iso[0] + ":" + t_mjd.iso[1])
    else:
        min_mjd = min(t["mjd"])
        max_mjd = max(t["mjd"])
        t_mjd = Time([min_mjd, max_mjd], scale = "utc", format = "mjd")
        mjd_info = "MJD Range: %0.2f:%0.2f\n" % (min_mjd, max_mjd)
        mjd_info += (t_mjd.iso[0] + ":" + t_mjd.iso[1])

    for band in bands:
        conv = convs[m]
        area_file = file[:-5] + "_area_mangle_out_" + band.upper()
        depth_file = file[:-5] + "_depth_mangle_out_" + band.upper()
        areas = []
        with open(area_file, "r") as af:
            for line in af:
                area = line.split()[0]
                if area[0] == "0":
                    areas.append(float(area))

        depths = []
        with open(depth_file, "r") as df:
            for line in df:
                depth = line.split()[0]
                if depth[0] <> "w":
                    depths.append(float(depth))

        areas = np.array(areas)*((180.0/np.pi)**2)
        depths = np.array(depths) + conv

        print("Total area:", np.sum(areas))

        bins = np.linspace(min(depths), max(depths), 101)
        areas_out = np.zeros(100)
        n = 0
        while n < len(bins)-1:
            ids = np.where( (depths >= bins[n]) & (depths < bins[n + 1]) )[0]
            areas_out[n] += np.sum(areas[ids])
            n += 1

        cum_areas = np.zeros(100)
        n = 0
        for area in areas_out[::-1]:
            cum_areas[n] = area + cum_areas[n-1]
            n += 1

        tot_area = np.sum(areas)
        half_area = tot_area/2.0
        n = 0
        while n < len(cum_areas)-1:
            if half_area > cum_areas[n] and half_area < cum_areas[n+1]:
                half_lim = bins[::-1][1:][n]
            n += 1
        print("Magnitude limit for half the area:", half_lim)

        half_lim = "%0.2f" % half_lim
        tot_area = "%0.2f" % tot_area

        if save_data:
            t = Table(data = [bins[::-1][1:], cum_areas/np.sum(areas), cum_areas], names = ["MAG_LIM_AB", "AREA_NORM", "AREA"])
            t.write("/data/sr525/VHS_area_depth_" + band + ".fits", overwrite = True)

        if norm:
            ax.plot(bins[::-1][1:], cum_areas/np.sum(areas), color = colours[m], label = band)
        else:
            ax.plot(bins[::-1][1:], cum_areas, color = colours[m], label = band)
        if band is not "Ks":
            info += (" " + band + ": " + tot_area + "  " + half_lim + "\n")
        else:
            info += (band + ": " + tot_area + "  " + half_lim)
        data.append([tot_area, half_lim])
        m += 1

    import matplotlib as mpl
    mpl.rc('text', usetex=True)
    info = r"\begin{tabular}{ccc} & Area & Mag. Lim. \\ & (Deg$^{2}$) & 50\% \\\hline Y & " + data[0][0] + "& " + data[0][1] + r"\\ J &" + data[1][0] + "&" + data[1][1] + r"\\ H & " + data[2][0] + "&" + data[2][1] + r"\\ K$_{s}$ &" + data[3][0] + "&" + data[3][1] + " \end{tabular}"
    #fig_table = ax.table(cellText=data, rowLabels=[r"J", r"H", r"K$_{s}$"], colLabels=["Area", "Mag. Lim. 50%"], bbox = [0.05, 0.02, 0.35, 0.25], colWidths=[0.12, 0.22], cellLoc = "center")
    #fig_table.auto_set_font_size(False)
    #fig_table.set_fontsize(12)
    #for key, cell in fig_table.get_celld().items():
    #    cell.set_linewidth(0)
    ax.text(0.03, 0.05, info, transform = ax.transAxes, bbox = dict(facecolor = "white", alpha = 0.7, boxstyle="square,pad=0.3"), fontsize = 14)
    #ax.text(0.01, 0.91, mjd_info, transform = ax.transAxes, bbox = dict(facecolor = "white", alpha = 0.7), family = "monospace")

    title_file = ""
    for letter in file:
        if letter <> "_":
            title_file += (letter)
        else:
            title_file += ("\\" + letter)
    print(file, title_file)
    plt.suptitle(title_file)
    plotid()
    if system == "Vega":
        ax.set_xlabel("Magnitude Limit [Vega]")
        ax.set_xlim(17, 22)
    elif system == "AB":
        ax.set_xlabel("Magnitude Limit [AB]", fontsize = 15)
        ax.set_xlim(18, 21.7)
    if norm:
        ax.set_ylabel("Normalised Cumulative Area")
    else:
        ax.set_ylabel("Cumulative Area / sq degrees", fontsize = 15)
    plt.legend()
    plt.show()


def tile_depth(file, system="AB",
               showplot=False,
               saveplot=True, plotdir=None):

    print('Reading:', file)
    t = Table.read(file)

    bands = ["Y", "J", "H", "Ks"]
    mag_lim_list = []
    info = "    Tiles Median\n"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colours = ["b", "g", "y", "r"]

    if file == "/data/vhs/vsa/dqc/VHSv20140517/vhs_vsa_dqc_tiles_fs_metadata.fits":
        #0.55 is 3 sigma to 5 sigma depth conversion
        if system == "AB":
            convs = [0.618-0.55, 0.937-0.55, 1.384-0.55, 1.839-0.55]
        else:
            convs = [-0.55, -0.55, -0.55, -0.55]
    else:
        if system == "AB":
            convs = [0.618, 0.937, 1.384, 1.839]
        else:
            convs = [0.0, 0.0, 0.0, 0.0]
    m = 0

    if "YMINRA" in t.columns:
        min_mjd = min(min(t["YMJDOBS"]), min(t["JMJDOBS"]), min(t["HMJDOBS"]), min(t["KSMJDOBS"]))
        max_mjd = max(max(t["YMJDOBS"]), max(t["JMJDOBS"]), max(t["HMJDOBS"]), max(t["KSMJDOBS"]))
        t_mjd = Time([min_mjd, max_mjd], scale = "utc", format = "mjd")
        mjd_info = "MJD Range: %0.2f:%0.2f\n" % (min_mjd, max_mjd)
        mjd_info += (t_mjd.iso[0] + ":" + t_mjd.iso[1])
    else:
        min_mjd = min(t["mjd"])
        max_mjd = max(t["mjd"])
        t_mjd = Time([min_mjd, max_mjd], scale = "utc", format = "mjd")
        mjd_info = "MJD Range: %0.2f:%0.2f\n" % (min_mjd, max_mjd)
        mjd_info += (t_mjd.iso[0] + ":" + t_mjd.iso[1])

    for band in bands:

        if band.upper() + "MINRA" in t.columns:
            ids = np.where( (t[band.upper() + "MINRA"] > -999999000.0) )[0]
            #Some of the data is at 3 sigma not 5
            mag_lim_list.append(t[band.upper() + "_DEPTH_DYE2006"][ids])

        else:
            ids = np.where( (t["filtname"] == band) )[0]
            print(len(ids))
            mag_lim_list.append(t["maglim"][ids])

    for (mags, band) in zip(mag_lim_list, bands):
        bins = np.linspace(min(mags), max(mags), 101)
        num_tiles = np.zeros(100)
        n = 0
        while n < len(bins)-1:
            ids = np.where( (mags >= bins[n]) & (mags < bins[n + 1]) )[0]
            num_tiles[n] += len(ids)
            n += 1

        cum_tiles = np.zeros(100)
        n = 0
        for num_tile in num_tiles[::-1]:
            cum_tiles[n] = num_tile + cum_tiles[n-1]
            n += 1

        plt.plot(bins[::-1][1:]+convs[m], cum_tiles/len(mags), label = band, color = colours[m])
        num_tiles = len(mags)
        med = np.median(mags)
        if band is not "Ks":
            info += (" " + band + ": " + str(num_tiles) + "  %0.2f\n" % (med + convs[m]))
        else:
            info += (band + ": " + str(num_tiles) + "  %0.2f" % (med + convs[m]))
        conv = convs[m]
        if conv <= 0.0:
            ax.set_xlabel("Magnitude Limit [Vega]")
        elif conv > 0.0:
            ax.set_xlabel("Magnitude Limit [AB]")
        if system == "AB":
            ax.set_xlim(18,22)
        else:
            ax.set_xlim(17, 22)

        m += 1
    ax.text(0.01, 0.02, info, transform = ax.transAxes, bbox = dict(facecolor = "white", alpha = 0.7), family = "monospace")
    ax.text(0.01, 0.90, mjd_info, transform = ax.transAxes, bbox = dict(facecolor = "white", alpha = 0.7), family = "monospace")
    plt.ylabel("Fraction with deeper than mag limit")
    plotid()
    plt.title(file, fontsize='medium')
    plt.legend(loc = "best")

    plotfile = 'tile_depth_hist_Vega.png'
    if system == 'AB':
        plotfile = 'tile_depth_hist_AB.png'

    if saveplot:
        print('Saving:', plotfile)
        plt.savefig(plotfile)

    if showplot:
        plt.show()

    return


def depth_compare():

    t = Table.read("/data/vhs/vistaqc_20140131_tiles_vhs.fits")
    t1 = Table.read("/data/vhs/vsa/dqc/VHSv20140517/vhs_vsa_dqc_tiles_fs_metadata.fits")

    for band in ["Y", "J", "H", "Ks"]:
        ids = np.where( (t["filtname"] == band) )[0]
        tb = t[ids]
        dists, inds = match_lists.match_lists(np.rad2deg(t1["RA"]), np.rad2deg(t1["DEC"]), tb["ra"], tb["dec"], 0.2)
        ids = np.where( (inds <> len(tb)) )[0]
        t11 = t1[ids]
        tb1 = tb[inds[ids]]
        plt.plot(tb1["maglim"], t11[band.upper() + "_DEPTH_DYE2006"], "k.")
        plt.xlabel("maglim")
        plt.ylabel("DEPTH_DYE2006")
        plt.title(band)
        plt.show()

def depth_check(file):

    t = Table.read(file)

    bands = ["Y", "J", "H", "Ks"]
    for band in bands:
        ids = np.where( (t["filtname"] == band) )[0]
        tb = t[ids]
        zpt = tb["magzpt"]- ((tb["amstart"]-1.0)*tb["extinct"])
        mag_check = zpt-2.5*np.log10(5.0*tb["skynoise"]*np.sqrt(1.0*np.pi*tb["rcore"]**2)/tb["exptime"])-tb["apcor3"]

        plt.plot(tb["maglim"], mag_check-tb["maglim"], "k.")
        plt.axhline(y = 0.0)
        plt.xlabel("CASU maglim")
        plt.ylabel("Calculated - CASU maglim")
        plt.title(file + " " + band)
        plt.show()

        print(np.median(mag_check-tb["maglim"]))

        plt.plot(tb["seeing"], mag_check-tb["maglim"], "k.")
        plt.xlabel("seeing")
        plt.ylabel("Calculated - CASU maglim")
        plt.title(file + " " + band)
        plt.show()


def getargs(verbose=False):
    """

    parse command line arguements

    not all args are active


    """
    import sys
    import pprint
    import argparse

    # there is probably a version function out there
    __version__ = '0.1'

    description = 'This is a template using getargs'
    epilog = """WARNING: Not all options may be supported
             """
    parser =  argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    # the destination defaults to the option parameter
    # defaul=False might not be needed

    parser.add_argument("--infile",
                        help="Input file name")

    parser.add_argument("--magsys",
                        default='AB',
                        help="Magnitude system [AB or Vega]")


    parser.add_argument("--configfile",
                        default=None,
                        help="configuration file")

    parser.add_argument("--debug",
                        action='store_true',
                        help="debug option")

    parser.add_argument("--verbose", default=verbose,
                        action='store_true',
                        help="verbose option")

    parser.add_argument("--version", action='store_true',
                        help="verbose option")


    args = parser.parse_args()


    if args.debug or args.verbose:
        print('Number of arguments:', len(sys.argv),
              'arguments: ', sys.argv[0])

    if args.debug or args.verbose:
        pprint.pprint(args)

    if args.version:
        print('version:', __version__)
        sys.exit(0)


    return args



if __name__ == '__main__':

    from librgm.getconfig import getconfig

    args = getargs()
    debug = args.debug

    configfile = 'vhs_mk_mangle_mask.cfg'
    if args.configfile is not None:
        configfile = args.configfile
    config = getconfig(configfile=configfile, debug=debug)

    inpath = config.get('CONFIG', 'inpath')
    infilename = config.get('CONFIG', 'infilename')
    infile = inpath + infilename
    print('infile:', infile)

    manglepath = config.get('CONFIG', 'manglepath')

    outpath = config.get('CONFIG', 'outpath')

    # file = "/data/vhs/vsa/dqc/VHSv20140517/vhs_vsa_dqc_tiles_fs_metadata.fits"

    # file = "/data/sr525/VHS/2016/vhs_vsa_dqc_tiles_fs_metadata.fits"
    #/data/sr525/VHS/2016 has various files vhs_vsa_dqc_tiles_fs_metadata...
    # infile = "/data/vhs/dqc/2014/vistaqc_20140131_tiles_vhs.fits"

    # plot cdf for tile depths
    tile_depth(infile, system = "AB", saveplot=True, plotdir='./')
    tile_depth(infile, system = "Vega", saveplot=True, plotdir='./')

    # this works
    area_file, depth_file = run_mangle(infile, "all", 0.0)

    print()
    print('area_file:', len(area_file))
    print('area_file:', area_file)
    print()
    print('depth_file:', len(depth_file))
    print('depth_file:', depth_file)
    print()

    print('run area_depth_all')
    print('infile:', infile)

    #area_depth_all(infile, system = "AB", norm = False)

    #sys.exit()

    #area_depth_all(infile, system = "AB", norm = True, save_data = True)

    #area_depth_all(file, system = "Vega", norm = False)
    #area_depth_all(file, system = "Vega", norm = True)

    #area_file, depth_file = run_mangle(infile, "all", 0.0)

    area_depth(infile, area_file, depth_file, "all", 0.0)


    """

    tile_depth(file, system = "AB")
    tile_depth(file, system = "Vega")
    area_depth_all(file, system = "AB", norm = False)
    area_depth_all(file, system = "AB", norm = True)
    area_depth_all(file, system = "Vega", norm = False)
    area_depth_all(file, system = "Vega", norm = True)

    #file = "/data/vhs/dqc/2014/vistaqc_20140131_tiles_vhs.fits"
    #area_depth_all(file)
    #tile_depth(file)
    #depth_compare()
    #depth_check(file)
    exptime = False
    #Need a factor of 2.5log10(3/5) = 0.55 to change 3 sigma limits to 5 sigma.
    #Takes it as Vega if conversion is < 0.6

    convs = [0.618-0.55, 0.937-0.55, 1.839-0.55, 1.384-0.55]
    #convs = [-0.55, -0.55, -0.55, -0.55]
    for (band, conv) in zip(["Y", "J", "H", "KS"], convs):
        area_file = file[:-5] + "_area_mangle_out_" + band
        depth_file = file[:-5] + "_depth_mangle_out_" + band
        area_file, depth_file = run_mangle(file, band, conv, exptime = exptime)
        area_depth(file, area_file, depth_file, band, conv)
    """
