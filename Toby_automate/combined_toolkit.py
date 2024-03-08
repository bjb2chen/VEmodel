import os
import shutil
import sys
import subprocess
import project_parameters as pp

def process_optimized_rhf(file):

    print("You have selected 1 - OPTIMIZED RHF, please ensure the file given to this script is GMS.dat")

    with open(file, 'r') as f:
        content = f.read()

    # Remove existing vec.dat
    if os.path.exists("vec.dat"):
        os.remove("vec.dat")

    # Extract section
    start = content.find("OPTIMIZED RHF")
    end = content.find("END", start)
    section = content[start:end+4] #include the $END, otherwise ends at $\n

    # Write to vec.dat
    with open("vec.dat", 'w') as f:
        f.write(section+"\n")

    print("OPTIMIZED RHF orbitals prepared in vec.dat file")

def process_optimized_rohf(file):

    print("You have selected 2 - OPTIMIZED ROHF, please ensure the file given to this script is GMS.dat")

    with open(file, 'r') as f:
        content = f.read()

    # Remove existing vec.dat
    if os.path.exists("vec.dat"):
        os.remove("vec.dat")

    # Extract section
    start = content.find("OPTIMIZED ROHF")
    end = content.find("END", start)
    section = content[start:end+4]

    # Write to vec.dat
    with open("vec.dat", 'w') as f:
        f.write(section+"\n")

    print("OPTIMIZED ROHF orbitals prepared in vec.dat file")

def process_semi_canonical(file):

    print("You have selected 3 - Semi-canonical MOs, please ensure the file given to this script is GMS.dat")

    with open(file, 'r') as f:
        content = f.read()

    # Remove existing active_space_orbs
    if os.path.exists("active_space_orbs"):
        os.remove("active_space_orbs")

    # Extract section
    start = content.find("Semi-canonical MOs")
    #print(start)
    end = content.find("END", start)
    #print(end)
    section = content[start:end+4]
    #print(section)

    # Write to active_space_orbs
    with open("active_space_orbs", 'w') as f:
        f.write(section+"\n")

    print("Semi-canonical MOs prepared in active_space_orbs file")

def process_mcscf_natural_orbitals(file):

    print("You have selected 4 - MCSCF Natural Orbitals, please ensure the file given to this script is GMS.dat")

    with open(file, 'r') as f:
        content = f.read()

    if os.path.exists("nat_orbs_mcscf"):
        os.remove("nat_orbs_mcscf")

    start = content.find("NATURAL ORBITALS OF MCSCF")
    end = content.find("END", start)
    section = content[start:end+4]

    with open("nat_orbs_mcscf", 'w') as f:
        f.write(section+"\n")

    print("Natural MCSCF orbitals prepared in nat_orbs_mcscf file")

def process_optimized_mcscf(file):

    print("You have selected 5 - OPTIMIZED MCSCF, please ensure the file given to this script is GMS.dat")

    with open(file, 'r') as f:
        content = f.read()

    if os.path.exists("optimized_mscsf"):
        os.remove("optimized_mscsf")

    start = content.find("OPTIMIZED MCSCF")
    end = content.find("END", start)
    section = content[start:end+4]

    with open("optimized_mscsf", 'w') as f:
        f.write(section+"\n")

    print("Natural MCSCF orbitals prepared in optimized_mscsf file")

def process_dmo_group(file):

    bash_script = f"""
    echo "You have selected 6 - DMO group, please ensure the file given to this script is GMS.dat"
    read -p $'Please input (x) number of DMOs and (y) the initial DMO in format: x y \nYour choice:' ndmo inidmo
    
    sed -n "/Semi-canonical/,/END/ p" {file} > ztemp
    
    norb=$(grep '^... 1' ztemp | wc -l)
    echo "number of orbitals: $norb"
    
    nline=$(wc -l ztemp | cut -d' ' -f1)
    nline_per_orb=$(echo "($nline - 3) / $norb" | bc -l | cut -d'.' -f1)
    echo "Number of lines for each orbita: $nline_per_orb"
    
    ini_line=$(echo "($inidmo - 1) * $nline_per_orb + 3" | bc -l)
    end_line=$(echo " $ini_line + $ndmo * $nline_per_orb - 1" | bc -l)
    
    echo $ini_line $end_line
    
    sed -n "$ini_line, $end_line p" ztemp > ztemp2
    
    # Remove existing dmo.dat
    if [ -f dmo.dat ]; then
      rm -f dmo.dat
    fi
    
    echo ' $DMO    ' >> dmo.dat
    for idmo in $(seq 1 1 $ndmo); do
      ini_line=$(echo "($idmo - 1) * $nline_per_orb" + 1 | bc -l)
      end_line=$(echo " $ini_line + $nline_per_orb - 1" | bc -l)
      if [ $idmo -lt 10 ]; then
        sed -n "$ini_line, $end_line p" ztemp2 | sed -e "s/^../ $idmo/g" >> dmo.dat
      else
        sed -n "$ini_line, $end_line p" ztemp2 | sed -e "s/^../$idmo/g" >> dmo.dat
      fi
    done
    
    echo ' $END' >> dmo.dat
    echo 'DMO group prepared in dmo.dat file'
    """
    os.system(bash_script)

def process_refdet_group(file):

    bash_script = f"""    
    echo "You have selected 7 - REFDET group, please ensure the file given to this script is GMS_DMOSTEP.out"
    nstate=`grep '# of states in CI      = ' {file}|tail -1|cut -d'=' -f2`   
    echo $nstate "states"
    
    sed -n "/DIABATIZATION FOR GMC-QDPT STATES/,/REFDET/ p" {file} > ztemp
    
    nstate_plus1=`echo "$nstate + 1"|bc -l`
    echo " State #    $nstate_plus1" >> ztemp
    
    
    # Remove the existing refdet.out
    if [ -f refdet.out ]
    then
      rm -f refdet.out
    fi
    
    # Remove the existing phase.out
    if [ -f phase.out ]
    then
      rm -f phase.out
    fi
    
    echo ' $REFDET     ' >> refdet.out
    echo " "$nstate >> refdet.out
    
    nleft=`expr $nstate % 5`
    nline_for_5=`echo "$nstate / 5"|bc -l`
    #echo $nline_for_5 $nleft
    
    for ist in $(seq 1 1 $nstate)
    do
     ist_plus1=`echo "$ist + 1"|bc -l`
     sed -n "/State..* $ist  Energy/,/State..* $ist_plus1  Energy/ p" ztemp |grep '^ .......[0-9]'|cut -c1-9,22-33 > ztemp2
     ndet=`wc -l ztemp2|cut -d' ' -f1`
     echo " $ist $ndet" >> refdet.out
     cat ztemp2 >> refdet.out
     igroup_5=0
     for jst in $(seq 1 1 $nstate)
     do
      modulo_jst_5=`expr $jst % 5`
    # echo $modulo_jst_5
      elemnt=' 0.00000000E+00'
      if [ $jst -eq $ist ]
      then
        elemnt=' 1.00000000E+00'
      fi
      if [ "$modulo_jst_5" -eq 1 ]
      then
        igroup_5=`echo " $igroup_5 + 1 "|bc -l`
        if [ $jst -ne $nstate ]
        then
          echo -n " $ist  $igroup_5$elemnt" >> phase.out
        else
          echo " $ist  $igroup_5$elemnt" >> phase.out
        fi
      elif [ "$modulo_jst_5" -eq 0 ] || [ $jst -eq $nstate ]
      then
        echo "$elemnt" >> phase.out
      else
        echo -n "$elemnt" >> phase.out
      fi
     done
    done
    
    cat phase.out >> refdet.out
    echo ' $END ' >> refdet.out
    echo "REFDET group prepared in refdet.out file"
    """

    os.system(bash_script)

def process_equilibrium_geometry(file):

    bash_script = f"""
    echo "You have selected 8 - EQUILIBRIUM GEOMETRY, please ensure the file given to this script is GMS_gh.out"
    # Remove existing ref_structure
    if [ -f ref_structure ]; then
      rm -f ref_structure
    fi
    # # Prompt user for the number of atoms
    # echo "Enter the number of atoms:"
    # read natoms
    
    # Extracting the table section
    table=$(sed -n "/EQUILIBRIUM GEOMETRY/,/INTERNUCLEAR DISTANCES/p" "{file}")
    
    # Extracting the last '{pp.Z}' lines from the table (excluding the last line of the section)
    echo "$table" | tail -n "$(({pp.Z} + 2))" | head -n -2  > "ref_structure"
    echo "EQUILIBRIUM GEOMETRY coordinates prepared in ref_structure"
    """

    os.system(bash_script)

if __name__ == "__main__":
    file = sys.argv[1]
    print(file)
    num = input("Do you want:\n (1) OPTIMIZED RHF \n (2) OPTIMIZED ROHF \n (3) Semi-canonical MOs \n (4) MCSCF Natural Orbitals \n (5) OPTIMIZED MCSCF \n (6) DMO group \n (7) REFDET GROUP \n (8) EQUILIBRIUM GEOMETRY \n Your choice: ")

    if num == "1":
        process_optimized_rhf(file)
    elif num == "2":
        process_optimized_rohf(file)
    elif num == "3":
        process_semi_canonical(file)
    elif num == "4":
        process_mcscf_natural_orbitals(file)
    elif num == "5":
        process_optimized_mcscf(file)
    elif num == "6":
        process_dmo_group(file)
    elif num == "7":
        process_refdet_group(file)
    elif num == "8":
        process_equilibrium_geometry(file)
    else:
        print("Invalid selection.")
