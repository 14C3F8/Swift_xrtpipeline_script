# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2  2023

@author: Songbo Gao
@Github page: https://github.com/14C3F8/Swift_xrtpipeline_script

"""

import os, glob, re, subprocess, time, shutil
import astropy.units as u
from astropy.coordinates import SkyCoord


def check_requirements():
    """
    Check if the source and back region files are in the current directory and if the HEASoft and CALDB environments are set up.
    Returns a tuple with the source region file and the back region file if all requirements are met, and (None, None) otherwise.
    """
    source_region = glob.glob('src.reg') # source region file
    back_region = glob.glob('back.reg') # background region file
    
    def is_headas_initialized():
        """
        Check if the HEASoft environment has been initialized by running the $HEADAS/headas-init.sh script.
        Returns True if the environment has been initialized, and False otherwise.
        """
        if 'HEADAS' not in os.environ:
            return False
        
        if not shutil.which('xrtversion'):
            return False
        
        return True
    
    if not source_region or not back_region:
        print("  \033[1;31mSource or back region files not found. Please check the current directory.\033[0m")
        return False, (None, None)

    if 'HEADAS' not in os.environ or 'CALDB' not in os.environ:
        print("  \033[1;31mPlease set up the HEASoft and CALDB environment variables.\033[0m")
        return False, (None, None)
    
    if not is_headas_initialized():
        print("  \033[1;31mPlease run the $HEADAS/headas-init.sh script to set up the HEASoft environment.\033[0m")
        return False, (None, None)
    
    return True, (source_region[0], back_region[0])

def confirm(prompt, default='yes'):
        """
        :param prompt: 
        :param default: yes
        """
        valid = {"yes": True, "y": True, "no": False, "n": False}
        if default is None:
            prompt += " [y/n] "
        elif default == "yes":
            prompt += " [Y/n] "
        elif default == "no":
            prompt += " [y/N] "
        else:
            raise ValueError(f"Invalid input.")

        while True:
            choice = input(prompt).lower()
            if default is not None and choice == '':
                return valid[default]
            elif choice in valid:
                return valid[choice]
            else:
                print("please input 'yes' or 'no' ('y' or 'n')")


def get_object_coordinates(object_name,region_file):
    
    libraries_installed = True
    try:
        from astroquery.simbad import Simbad
        import concurrent
        from concurrent.futures import ThreadPoolExecutor
    except ImportError:
        libraries_installed = False
    
    def query_object_coordinates(object_name):
        result_table = Simbad.query_object(object_name)
        ra = result_table['RA'].data[0]
        dec = result_table['DEC'].data[0]
        coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        return coord.ra.deg, coord.dec.deg

    def parse_reg_file(region_file):
        with open(region_file, 'r') as file:
            lines = file.readlines()

        line = lines[3]
        pattern = r'\((.+?)\,(.+?)\,.*\)'
        match = re.search(pattern, line)
        ra_str, dec_str = match.group(1).strip(), match.group(2).strip()

        coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
        ra_deg, dec_deg = coord.ra.deg, coord.dec.deg

        return ra_deg, dec_deg    
    
    if libraries_installed:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future = executor.submit(query_object_coordinates, object_name)
            try:
                print('  Get source position from SIMBAD...\n')
                ra_deg, dec_deg = future.result(timeout=20)
            except concurrent.futures.TimeoutError:
                print('  \033[1;30mTime out\033[0m, reading source position from region file.\n')
                ra_deg, dec_deg = parse_reg_file(region_file)
    else:
        print('  \033[1;30mAstroquery or futures package not installed\033[0m, reading source postion from region file.\n')
        ra_deg, dec_deg = parse_reg_file(region_file)

    return ra_deg, dec_deg 


def get_observation_ids(path):
    """
    Get the observation IDs from the /path/to/data directory.
    Returns a list of observation IDs.
    """
    
    observation_dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    observation_ids = [os.path.basename(d) for d in observation_dirs]
    return observation_ids, len(observation_dirs)


def check_and_create_result_directory(path):
    """
    Check if the given directory exists. If not, create it.
    """
    if not os.path.exists(path):
        print('  result directory not exists, creating a new directory instead.')
        os.makedirs(path)

def find_event_files(obsid, result_dir):
    """
    find event files generated from xrtpipeline.
    """
    pc_event_file_pattern = os.path.join(result_dir, obsid, f'sw{obsid}xpc*po_cl.evt')
    wt_event_file_pattern = os.path.join(result_dir, obsid, f'sw{obsid}xwt*po_cl.evt')

    pc_event_files = glob.glob(pc_event_file_pattern)
    wt_event_files = glob.glob(wt_event_file_pattern)

    if not pc_event_files:
        pc_event_files = None

    if not wt_event_files:
        wt_event_files = None

    return pc_event_files, wt_event_files   

def run_xrtpipeline(obsid,ra,dec,data_path,result_path):
    """
    run xrtpipeline
    """
    error_flag = True
    start_time = time.time()
    log_filename = f"./log/{obsid}_step1_log.txt"
    error_filename = "step1_error_obsids.txt"

    with open(log_filename, "w") as log_file:
        # Run xrtpipeline
        
        xrtpipeline_cmd = f"punlearn xrtpipeline && xrtpipeline indir={data_path}/{obsid} outdir={result_path}/{obsid} steminputs=sw{obsid} exitstage=2 srcra={ra} srcdec={dec} createexpomap=yes cleanup=no clobber=yes"
        env_vars = os.environ.copy()
        print(f'  running xrtpipeline for {obsid}...')
        subprocess.run(xrtpipeline_cmd, shell=True, env=env_vars, stdout=log_file, stderr=log_file)
        
    pc_event_file, wt_event_file = find_event_files(obsid, result_path)
        
    if pc_event_file is None:
        print('  No PC mode pointing data found!')
        
    if wt_event_file is None:
        print('  No WT mode pointing data found!')
            
    if pc_event_file is None and wt_event_file is None:
        with open(error_filename, "a") as error_file:
            error_file.write(f"{obsid}_pc\n")
        print(f"  \033[1;31mError when processing observation {obsid} data with xrtpipeline. Please check your data or logfile.\033[0m")
        error_flag = False
 
    time_elapsed = time.time() - start_time
    print(f"  Finished running xrtpipeline for {obsid} in {time_elapsed:.2f} seconds.")

    return error_flag

def find_xrtproducts_files(obsid, result_dir):
    """
    Find event files generated from xrtpipeline.
    """
    pc_spectrums = []
    wt_spectrums = []
    pc_lightcurves = []
    wt_lightcurves = []
    
    obsid_result_dir = os.path.join(result_dir, obsid)
    
    pc_spectrum_pattern = os.path.join(obsid_result_dir, f'sw{obsid}xpc*posr.pha')
    wt_spectrum_pattern = os.path.join(obsid_result_dir, f'sw{obsid}xwt*posr.pha')
    pc_lightcurve_pattern = os.path.join(obsid_result_dir, f'sw{obsid}xpc*po_corr.lc')
    wt_lightcurve_pattern = os.path.join(obsid_result_dir, f'sw{obsid}xwt*po_corr.lc')
    
    pc_spectrums = glob.glob(pc_spectrum_pattern)
    wt_spectrums = glob.glob(wt_spectrum_pattern)
    pc_lightcurves = glob.glob(pc_lightcurve_pattern)
    wt_lightcurves = glob.glob(wt_lightcurve_pattern)
    
    if not pc_spectrums:
        pc_spectrums = None
    
    if not wt_spectrums:
        wt_spectrums = None
            
    if not pc_lightcurves:
        pc_lightcurves = None
    
    if not wt_lightcurves:
        wt_lightcurves = None        
            
    return pc_spectrums, wt_spectrums, pc_lightcurves, wt_lightcurves   

def run_xrtproducts(obsid,data_path,result_path,event_file,src_reg,bkg_reg,binsize,mode,pilow,pihigh):
    """
    Run xrtproducts
    """
    error_flag = True
    start_time = time.time()
    log_filename = f"./log/{obsid}_{mode}_step2_log.txt"
    error_filename = "step2_error_obsids.txt"

    with open(log_filename, "w") as log_file:
        xrtproducts_cmd = f"xrtproducts infile={event_file} regionfile={src_reg} bkgextract=yes bkgregionfile={bkg_reg} outdir={result_path}/{obsid} stemout=DEFAULT binsize={binsize} pilow={pilow} pihigh={pihigh} expofile={data_path}/{obsid}/xrt/products/sw{obsid}x{mode}_ex.img.gz attfile={data_path}/{obsid}/auxil/sw{obsid}pat.fits.gz hdfile={data_path}/{obsid}/xrt/hk/sw{obsid}xhd.hk.gz"
        env_vars = os.environ.copy()
        print(f'  running xrtproudcts for {obsid}...')
        subprocess.run(xrtproducts_cmd, shell=True, env=env_vars, stdout=log_file, stderr=log_file)
        
    pc_spectrums, wt_spectrums, pc_lightcurves, wt_lightcurves = find_xrtproducts_files(obsid, result_path)
    if pc_spectrums is None and wt_spectrums is None and pc_lightcurves is None and wt_lightcurves is None:
        with open(error_filename, "a") as error_file:
            error_file.write(f"{obsid}_pc\n")
        print(f"  \033[1;31mError when processing observation {obsid} data with xrtproducts. Please check your data or logfile.\033[0m")
        error_flag = False
    
    
    time_elapsed = time.time() - start_time
    print(f"  Finished running xrtproducts for {obsid} {mode} data in {time_elapsed:.2f} seconds.")
    
    return error_flag

def find_response(mode,pcrmf,wtrmf):
      
    caldb_path = os.getenv('CALDB')
    rmf_path = os.path.join(caldb_path, "data/swift/xrt/cpf/rmf")
    
    if mode == 'pc':
        rmf = os.path.join(rmf_path,pcrmf)
    elif mode == 'wt':
        rmf = os.path.join(rmf_path,wtrmf)
    
    return rmf

def run_xrtmkarf(obsid,data_path,result_path,mode,spectrum_file,pcrmf,wtrmf):
    """
    Run xrtmkarf
    """
    
    error_flag=True
    start_time = time.time()
    log_filename = f"./log/{obsid}_{mode}_step3_log.txt"
    error_filename = "step2_error_obsids.txt"
    
    resp_file = find_response(mode,pcrmf,wtrmf)
    spec_pattern = spectrum_file[-10:-4]
    
    with open(log_filename, "w") as log_file:
        xrtmkarf_cmd=f"xrtmkarf outfile={result_path}/{obsid}/{obsid}_{mode}{spec_pattern}.arf phafile={spectrum_file} expofile={data_path}/{obsid}/xrt/products/sw{obsid}x{mode}_ex.img.gz rmffile={resp_file} srcx=0 srcy=0 psfflag=yes"
        env_vars = os.environ.copy()
        print(f'  running xrtmkarf for {obsid}...')
        subprocess.run(xrtmkarf_cmd, shell=True, env=env_vars, stdout=log_file, stderr=log_file)
        
    arf_file = f'{result_path}/{obsid}/{obsid}_{mode}{spec_pattern}.arf'
    arf_pattern = glob.glob(arf_file)
    if not arf_file:
        with open(error_filename, "a") as error_file:
            error_file.write(f"{obsid}_pc\n")
        print(f"  \033[1;31mFailed generate arf for {obsid} {mode} spectrum using xrtmkarf. Please check your data or logfile.\033[0m")
        error_flag = False
        arf_file = None
    
    time_elapsed = time.time() - start_time
    print(f"  Finished running xrtmkarf for {obsid} {mode} spectrum {spectrum_file[-24:]} in {time_elapsed:.2f} seconds.")
    
    return error_flag, arf_pattern, resp_file

def run_grppha(obsid,result_path,spectrum_binsize,src_spectrum_file,resp_file,arf_file):
    
    spec_pattern = src_spectrum_file[-24:-4]
    back_spectrum = src_spectrum_file[:-6]+'bkg.pha'
    arf_file = arf_file[0]
    grppha_cmd = f"grppha infile={src_spectrum_file} outfile={result_path}/{obsid}/{spec_pattern}_grp.pha chatter=0 comm='group min {spectrum_binsize} & chkey RESPFILE {resp_file} & chkey ANCRFILE {arf_file} & chkey BACKFILE {back_spectrum} & exit'"
    env_vars = os.environ.copy()
    subprocess.run(grppha_cmd,shell=True, env=env_vars,stdout=subprocess.DEVNULL)
    print(f'  Finished grouping PHA files. Output file is {spec_pattern}_grp.pha.')
