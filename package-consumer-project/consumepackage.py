#from __future__ import print_function
import os,sys
import numpy as np

#from pytec import stations as st
#from pytec import tec

f_path = ".env"
pwd = os.environ["PWD"]

if os.path.exists(f_path):
    f = open(f_path,'r')
    PYTEC_PATH = f.readlines()[0]
    print ("PATH used to stored GNSS and receiver bias already set to: "+ PYTEC_PATH)
    f.close()
else:
    PYTEC_PATH = pwd + "/.pytec"
    print ("Please define PYTEC directory where all intermediate data will be stored, this can be an external disk (consider at least 1GB)")
    print ("    This folder will contain two sub folders:")
    print ("       -- TEC: containing the bias of the receptors")
    print ("       -- GPS: saving the position of each satellite to avoid reprocessing")
    ans = str(input ("Use default directory '"+PYTEC_PATH+"'? (y/n) "))

    while True:
        if ans == "y":
            f = open(f_path,"w")
            home = os.environ["HOME"]
            f.write("PYTEC_PATH=\""+PYTEC_PATH+"\"")
            f.close()
            print ("pytec directory set up to "+PYTEC_PATH)
            if not os.path.exists(PYTEC_PATH):
                try: os.mkdir(PYTEC_PATH)
                except OSError as e:
                    if e.errno!=17: 
                        print ("FAIL creation of directory "+PYTEC_PATH, e )
                        print ("Installation failed, please check that  folder ",PYTEC_PATH, "can be created")
            break
        elif ans == "n":
            PYTEC_PATH = str(input("Please enter the path of the directory you want to use (it must exist, it will NOT be created)\n"))
            if not os.path.exists(PYTEC_PATH):
                print (PYTEC_PATH,"does not exist, restart installation with an existing folder (or use default).\n")
                print ("ERROR: Installation failed")
                break
            f = open(f_path,"w")
            home = os.environ["HOME"]
            f.write("PYTEC_PATH=\""+PYTEC_PATH+"\"")
            f.close()
            print ("pytec directory set up to "+PYTEC_PATH)
            break
        else:
            print ("Did not understand, please answer 'y' or 'n'")
            ans = str(input ("Use default directory '"+PYTEC_PATH+"'? (y/n) "))
            
        
    print ("You can change it by editing or removing the file "+f_path+", and reinstalling pytec")

sys.exit()

#if not os.environ.get("PYTEC_PATH"):
if "PYTEC_PATH" not in os.environ:
    print ("First of all you need to define the environment variable 'PYTEC_PATH' pointing to the directory where all intermediate data will be stored, can be an external disk (consider at least 1GB)")
    print ("    This folder will contain two sub folders:")
    print ("       -- TEC: containing the bias of the receptors")
    print ("       -- GPS: saving the position of each satellite to avoid reprocessing")
    print ("\n Note that if this environment variable is not defined pytec can work fine, but will not avoid reprocessing and computing TEC could take more time")
    ans = str(input ("Do you want to define it? (y/n) (Recomended) "))
    
    if ans=='y':
        home = os.environ["HOME"]
        while True:
            path = input("Which path for pytec? (Press Enter to use default "+ home+"/.pytec)" )
            if path=='':
                path = home+"/.pytec"
            if not os.path.exists(path):
                try: os.mkdir(path)
                except OSError as e:
                    if e.errno!=17: 
                        print ("FAIL creation of directory "+path, e )
                        again = str(input ("keep trying? y/n"))
                        if again == 'n': 
                            print ("No fixed path will be defined for pytec this can be done by hand later in the ~/.bashrc file)")
                            break
                else: 
                    print ("Successfully created the directory "+path)
                    bashrc = open(home +'/.bashrc', 'a')
                    bashrc.write("\nexport PYTEC_PATH="+path)
                    bashrc.close()
                    break
            else:
                print (path + " exists, will be used for pytec use")
                bashrc = open(os.environ.get("HOME")+'/.bashrc', 'a')
                bashrc.write("\nexport PYTEC_PATH="+path)
                bashrc.close()
                break
                
    else:
        print ("No fixed path will be defined for pytec this can be done by hand later in the ~/.bashrc file)")
else:
    pytec_path = os.environ.get("PYTEC_PATH")
    print (pytec_path + "  exists")
    if not os.path.exists(pytec_path):
        try: os.mkdir(pytec_path)
        except OSError as e:
            if e.errno!=17: 
                print ("FAIL creation of directory "+pytec_path, e )
                print ("This is inconsistent with the fact that the PYTEC_PATH environment variable exists, but you can set this up:")
                print ("    1 -- Remove the PYTEC_PATH variable from the ~/.bashrc file")
                print ("    2 -- Run 'unset PYTEC_PATH'")
                print ("    3 -- Run 'source ~/.bashrc'")
                print ("    4 -- Rerun this script")

        else: 
            print ("Successfully created directory "+pytec_path)

        
