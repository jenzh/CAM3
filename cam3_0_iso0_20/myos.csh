#!/bin/csh -f
#
# Lists files which have been modified from the original.
#   Files given with no path and with a .o extension, so that they may be removed to 
#   add debugging (e.g.) options in the "bld" directorty
#
set verb = 0		# verbose outout
set nl = '-n' 
if ($verb > 0) set nl = ' ' 


#set srcori = '/data2/csm/cam3_0/cam1/models'		# originaal code
set srcori = '../cam3_0_iso0_11/cam1/models'		# older isotope version
set srciso = 'cam1/models'				# this source tree

#
#
#
set diffops = '-a -b -B -i'

#
# This one all diffs
#
diff $diffops -r $srcori $srciso

#
# Do minimal outpout
#
if ($verb > 0) diff $diffops -r --brief $srciso $srcori
echo

#
# Grab the new files
#
set files = `diff $diffops -r --brief $srciso $srcori | grep F90 | grep ^'Only in' | grep -v $srcori | cut -f 2 -d':'`

if ($verb > 0) echo NEW FILES
foreach file ($files)
  set this = `echo $file | sed s/\.F90/\.o/g`
  echo $nl "$this "
end

#
# And the modified files
#
set files = `diff $diffops -r --brief $srciso $srcori | grep F90 | grep differ | cut -f 2 -d' '`

if ($verb > 0) echo MODIFIED FILES
foreach file ($files)
  set this = `basename $file | sed s/\.F90/\.o/g`
  echo $nl "$this "
end

echo
echo

# Done



