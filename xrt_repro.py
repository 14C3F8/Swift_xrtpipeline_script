# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 01:47:28 2023
@author: Songbo Gao
@Github page: https://github.com/14C3F8/14c3f8-astro-tools/Swift_XRT

"""

import os, sys, shutil, time, datetime, argparse
from xrt_repro_func import check_requirements, confirm, check_and_create_result_directory
from xrt_repro_func import get_observation_ids, find_event_files,run_xrtpipeline, find_xrtproducts_files
from xrt_repro_func import run_xrtproducts, run_xrtmkarf, run_grppha


def main(indir,outdir,obsobject,pc_binsize,wt_binsize,pilow,pihigh,pcrmf,wtrmf,spec_binsize):
    """
    Main function to run the script.
    """ 
    print()
    print(' SWIFT XRT DATA REDUCTION SCRIPT '.center(60, '-'))
    print()
    print(' Author: Songbo Gao '.center(60, ' '))
    print()
    print('Shandong Key Laboratory of Optical Astronomy'.center(60, ' ')) 
    print('and Solar-Terrestrial Environment,'.center(60, ' ')) 
    print('School of Space Science and Physics,' .center(60, ' '))
    print('Institute of Space Sciences, Shandong University'.center(60, ' '))
    print()
    print(' Github page: https://github.com/14C3F8/14c3f8-astro-tools/Swift_XRT '.center(60, ' '))
    print(' Email: 202117741@mail.sdu.edu.cn '.center(60, ' '))
    print()
    print(' Script created on April 2023'.center(60, ' '))
    print()
    print('--'.center(60, '-'))
    time.sleep(2)
    print(' Checking requirements '.center(60, '-'))
    print()
    
    if_requirements_satisfied, (src, back) = check_requirements()
    
    if if_requirements_satisfied:
        print('  All requirements sactified.\n')
        print(f'  Find region file {src} and {back}.\n')
    else:
        exit(1)
     
    result_directory = outdir
    check_and_create_result_directory(result_directory)
    
    log_path = os.path.join(outdir, 'log')
    if os.path.exists(log_path):
        shutil.rmtree(log_path)
        os.makedirs(log_path)
        print('  Logfile path exits, previous log files are deleted.\n')
    else:
        os.makedirs(log_path)
        print(f"  Logfile directory {log_path} created.\n")

    data_directory = indir
    observation_ids, total_observations = get_observation_ids(data_directory)
       
    from xrt_repro_func import get_object_coordinates
    
    ra, dec = get_object_coordinates(obsobject,src)
    src_ra = f"{ra:.3f}"
    src_dec = f"{dec:.3f}"
    time.sleep(2)
    dt = datetime.datetime.now()
    print(' INFORMATION '.center(60, '-'))
    print(f'''    
    Current time: \033[1;30m{dt.strftime('%Y-%m-%d %H:%M:%S')}\033[0m.
    
    Get object name \033[1;30m{obsobject}\033[0m.
    Find {total_observations} observations in {data_directory}.    
    Source position: \033[1;30mRA={src_ra} DEC={src_dec}\033[0m.
    Using region files: \033[1;30m{src} & {back}\033[0m to extract lightcurve & spectrum.
    
    \033[1;30mPhoton Counting (PC)\033[0m mode:
    - Using time binsize \033[1;30m{pc_binsize}s\033[0m.
    - Response matrix file used for mkarf: \033[1;30m{pcrmf}\033[0m.
    \033[1;30mWindowed Timing (WT)\033[0m mode:
    - Using time binsize \033[1;30m{wt_binsize}s\033[0m.
    - Response matrix file used for mkarf: \033[1;30m{wtrmf}\033[0m.
    
    Spectrum will be groupped using grppha with minium binsize \033[1;30m{spec_binsize}\033[0m.    
    ''')
    print('--'.center(60, '-'))
    print()
    
    if not confirm('  Are these parameters right?'):
        exit(1)    
    
    print()
    print('--'.center(60, '-'))
    print(f'''
    Ready to reprocess XRT data, please wait...
    Log files located at \033[1;30m{log_path}\033[0m. 
    \033[0;31mPlease check if there is any problem after running xrt_repro.py.\033[0m
    
    ''')
    print(' DATA REDUCTION '.center(60, '-'))
    print()
    time.sleep(1)
    error_obsid_file = 'error_obsid.txt'
    processed_count = 0
    error_obsid_count = 0

    for obsid in observation_ids:
        
        time.sleep(1)
        print(f'  {obsid}  '.center(60, '-'))
        print()
        Flag = True
        pc_step_2 = True
        wt_step_2 = True
        pc_step_3 = True
        wt_step_3 = True

        step_1=run_xrtpipeline(obsid,src_ra,src_dec,data_directory,result_directory)
        
        pc_event_files, wt_event_files = find_event_files(obsid,result_directory)
        
        if step_1:
            for pcevtfile in pc_event_files or []:
                pc_step_2 = run_xrtproducts(obsid,data_directory,result_directory,pcevtfile,src,back,pc_binsize,'pc',pilow,pihigh)
                if not pc_step_2:
                    Flag = False
        
            for wtevtfile in wt_event_files or []:
                wt_step_2 = run_xrtproducts(obsid,data_directory,result_directory,wtevtfile,src,back,wt_binsize,'wt',pilow,pihigh)
                if not wt_step_2:
                    Flag = False
        
        else:
            Flag = False
               
        pc_spectrums, wt_spectrums, pc_lightcurves, wt_lightcurves = find_xrtproducts_files(obsid, result_directory)
                  
        for pcspec in pc_spectrums or []:
            pc_step_3, pc_arf_name, pc_rmf_name = run_xrtmkarf(obsid,data_directory,result_directory,'pc',pcspec,pcrmf,wtrmf)
            if pc_step_3:
                run_grppha(obsid,result_directory,spec_binsize,pcspec,pc_rmf_name,pc_arf_name)
            else:
                Flag = False
             
        for wtspec in wt_spectrums or []:
            wt_step_3, wt_arf_name, wt_rmf_name = run_xrtmkarf(obsid,data_directory,result_directory,'wt',wtspec,pcrmf,wtrmf)
            if wt_step_3:
                run_grppha(obsid,result_directory,spec_binsize,wtspec,wt_rmf_name,wt_arf_name)
            else:
                Flag = False
                    
        if Flag is False: 
            error_obsid_count += 1
            with open(error_obsid_file, "a") as error_obsid:
                error_obsid.write(f"{obsid}\n")
        
        processed_count += 1
        print()
        print('  Finished  '.center(60, '-'))
        print(f"\n  Reprocessed \033[1;30m{processed_count}/{total_observations}\033[0m observations. \033[1;30m{total_observations - processed_count}\033[0m left.\n")
                
    print('   REPORT   '.center(60, '-'))
    print()
    print(f'  Finished running xrt_repro.py.'.center(60,' ' ))
    dt = datetime.datetime.now()
    print(dt.strftime('%Y-%m-%d %H:%M:%S').center(60,' '))
    print()
    print(f'  \033[1;30m{processed_count}\033[0m success, \033[1;30m{error_obsid_count}\033[0m failed.')
    
    print('  \033[1;31mDo not forget to check the log files.\033[0m')
    if error_obsid_count != 0:
        print('  \033[1;31mPlease check error_obsid.txt for the error observations.\033[0m')
    print()
    print('--'.center(60, '-'))

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Script for processing SWIFT XRT data.")    
    parser.add_argument("--indir", required=True, help="Path to the data directory.")
    parser.add_argument("--outdir", required=True, help="Path to the result directory.")
    parser.add_argument("--obsobject", required=True, help="Name of the source.")
    parser.add_argument("--pc-binsize", required=False, default=10, help="Time bin size for PC mode (s).The default value is 10s.")
    parser.add_argument("--wt-binsize", required=False, default=1, help="Time bin size for WT mode (s). The default value is 1s.")
    parser.add_argument("--pilow", required=False, default=20, help="Minium PI for lc extract. The default value is 20.")
    parser.add_argument("--pihigh", required=False, default=1000, help="Maxium PI for lc extract. The default value is 1000.")
    parser.add_argument("--pcrmf", required=False, default="swxpc0to12s6_20130101v014.rmf", help="PC response matrix. The default value is[swxpc0to12s6_20130101v014.rmf].")
    parser.add_argument("--wtrmf", required=False, default="swxwt0to2s6_20131212v015.rmf", help="WT response matrix. The default value is[swxwt0to2s6_20131212v015.rmf].")
    parser.add_argument("--spec_binsize", required=False, default=1, help="Minium spectrum photon number per bin. The default value is 1.")
    args = parser.parse_args()

    main(args.indir, args.outdir, args.obsobject, args.pc_binsize, args.wt_binsize, args.pilow, args.pihigh, args.pcrmf, args.wtrmf, args.spec_binsize)
