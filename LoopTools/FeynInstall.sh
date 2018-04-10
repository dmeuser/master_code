#! /bin/sh
# installation script for FeynArts, FormCalc, and LoopTools
# last modified 10 Feb 18 th

# if started from within a directory DIR (e.g. the user's
# home directory), the final directory structure will be:
#	DIR/$fadir
#	DIR/FeynArts (link to $fadir)
#	DIR/$fcdir
#	DIR/FormCalc (link to $fcdir)
#	DIR/$ltdir
#	DIR/LoopTools (link to $ltdir)

fadir=FeynArts-3.9
fatar=$fadir.tar.gz

fcdir=FormCalc-9.6
fctar=$fcdir.tar.gz

ltdir=LoopTools-2.13
lttar=$ltdir.tar.gz

cwd="`pwd`"

id=`id -u`

system="`uname -s`"

unset http_proxy

dltool() {
  wget="curl --location --remote-name --user-agent FeynInstall"
  $wget --version > /dev/null 2>&1 && return
  wget="wget --user-agent FeynInstall"
  $wget --version > /dev/null 2>&1 && return
  echo "Please install either curl or wget"
  exit 1
}

dltool


echo "Install FeynArts in $cwd/FeynArts?"
read yesno
case "$yesno" in
[yY]*)
  ifa=1
  make -f - << _EOF_ > install.log-FeynArts
InstallFeynArts:
	@echo "... downloading $fatar" 1>&2
	$wget http://feynarts.de/$fatar 2>&1
	rm -fr $fadir
	@echo "... unpacking tar file" 1>&2
	gunzip -c $fatar | tar xf -
	rm -f FeynArts
	ln -s $fadir FeynArts
	rm -f $fatar
	@echo "... done" 1>&2
_EOF_
  if test $? -eq 0 ; then
    echo "FeynArts installed successfully"
  else
    echo "Installation error, please check install.log-FeynArts"
  fi
  ;;
esac
echo ""


echo "Install LoopTools in $cwd/LoopTools?"
read yesno
case "$yesno" in
[yY]*)
  ilt=1
  make -f - << _EOF_ > install.log-LoopTools
InstallLoopTools:
	@echo "... downloading $lttar" 1>&2
	$wget http://www.feynarts.de/looptools/$lttar 2>&1
	rm -fr $ltdir
	@echo "... unpacking tar file" 1>&2
	gunzip -c $lttar | tar xf -
	@echo "... compiling" 1>&2
	(cd $ltdir && ./configure && \$(MAKE) default install clean) 2>&1
	@cd $ltdir && case "\`grep '^FC =' makefile\`,$system" in \\
	  *ifort*)	mf=makefile.quad-ifort ;; \\
	  *xlf*)	mf=makefile.quad-xlf   ;; \\
	  *f77*OSF1)	mf=makefile.quad-alpha ;; \\
	  *)		exit 0 ;; \\
	esac && \\
	{ echo "... compiling quadruple-precision version" 1>&2 ; \\
	  \$(MAKE) -f \$\$mf install clean 2>&1 ; }
	rm -f LoopTools
	ln -s $ltdir LoopTools
	rm -f $lttar
	@echo "... done" 1>&2
_EOF_
  if test $? -eq 0 ; then
    echo "LoopTools installed successfully"
  else
    echo "Installation error, please check install.log-LoopTools"
  fi
  ;;
esac
echo ""


echo "Install FormCalc in $cwd/FormCalc?"
read yesno
case "$yesno" in
[yY]*)
  ifc=1
  make -f - << _EOF_ > install.log-FormCalc
InstallFormCalc:
	@echo "... downloading $fctar" 1>&2
	$wget http://www.feynarts.de/formcalc/$fctar 2>&1
	rm -fr $fcdir
	@echo "... unpacking tar file" 1>&2
	gunzip -c $fctar | tar xf -
	rm -f FormCalc
	ln -s $fcdir FormCalc
	@echo "... compiling" 1>&2
	cd $fcdir && ./compile 2>&1
	rm -f $fctar
	@echo "... done" 1>&2
_EOF_
  if test $? -eq 0 ; then
    echo "FormCalc installed successfully"
  else
    echo "Installation error, please check install.log-FormCalc"
  fi
  ;;
esac
echo ""


mathcmd=math
shopt -s nullglob > /dev/null 2>&1
set --
case "$system" in
Darwin)
	mathcmd=MathKernel
	set -- /Applications/Mathematica*/Contents/MacOS \
	       $HOME/Desktop/Mathematica*/Contents/MacOS ;;
CYG*)
	set -- "$PROGRAMFILES/Wolfram Research/Mathematica"/* ;;
esac
for dir in "$@" ; do
  path="$dir:$path"
done
mathcmd="`PATH=\"$PATH:$path\" which $mathcmd`"

if "$mathcmd" -run "Print[7 673]; Exit" < /dev/null | grep 4711 > /dev/null ; then
  eval -- `"$mathcmd" -run '
    org[$Failed] = "";
    org[file_] := File /. FileInformation[file];
    Print["fa=\"" <> org[System\`Private\`FindFile["FeynArts\`"]] <> "\""];
    Print["fc=\"" <> org[System\`Private\`FindFile["FormCalc\`"]] <> "\""];
    Print["lt=\"" <> org[System\`Private\`FindFile["LoopTools"]] <> "\""];
    Exit[]
  ' < /dev/null | tail -2 | tr '\r' ' '`

  ask() {
    :
  }

  ask1() {
    dir="$cwd/$5"
    test "$dir/$3" -ef "$4" && return
    echo "Do you want to add $dir to Mathematica's \$Path,"
    echo "such that $1 can be loaded with just '$2'?"
    read yesno
    case "$yesno" in
    [yY]*)
	case "$system" in
	CYG*)	dir=`cygpath -w "$dir"`
		mmapath="$mmapath, \"${dir//\\/\\\\}\"" ;;
	*)	mmapath="$mmapath, \"$dir\"" ;;
	esac ;;
    esac
    echo ""
  }

  ask$ifa FeynArts "<< FeynArts\`" FeynArts.m "$fa" FeynArts
  ask$ifc FormCalc "<< FormCalc\`" FormCalc.m "$fc" FormCalc
  ask$ilt LoopTools "Install[\"LoopTools\"]" LoopTools "$lt" LoopTools/*/bin

  test -n "$mmapath" && "$mathcmd" -run "mmapath={0$mmapath}" -run '
    prefdir = ToFileName[$PreferencesDirectory, "Kernel"];
    If[ FileType[prefdir] === None, CreateDirectory[prefdir] ];
    hh = OpenAppend[ToFileName[prefdir, "init.m"]];
    Block[ {home = ToFileName[$HomeDirectory], $HomeDirectory, ToFileName},
      striphome[dir_] := If[ # == dir, #, ToFileName[$HomeDirectory, #] ]& @
        StringReplace[dir, home -> ""];
      SetAttributes[Write, HoldRest];
      With[ {path = striphome/@ Rest[mmapath]},
        Write[hh, $Path = Join[path, $Path]] ]
    ];
    Print["Modified ", Close[hh]];
    Exit[]
  ' < /dev/null | tail -1
else
  echo "Cannot run Mathematica (license problems?)."
  echo "Skipping modification of \$Path."
fi


cat << \_EOF_

-------------------------------------------------------------------

Thank you for using FeynArts, FormCalc, and LoopTools.

If you find any bugs, or want to make suggestions, or just write
fan mail, address it to Thomas Hahn <hahn@feynarts.de>.

Considering the manpower that has gone and still goes into the 
development of these packages, it is about fair that you cite the
following references if you use FeynArts, FormCalc, or LoopTools
to produce published results:

FeynArts 3:
  T. Hahn, Comput. Phys. Commun. 140 (2001) 418
  [hep-ph/0012260]

The MSSM model file of FeynArts:
  T. Hahn, C. Schappacher, Comput. Phys. Commun. 143 (2002) 54
  [hep-ph/0105349]
including counter-terms:
  T. Fritzsche, T. Hahn, S. Heinemeyer, H. Rzehak, C. Schappacher,
  Comput. Phys. Commun. 185 (2014) 1529 [arXiv:1309.1692]

FormCalc and LoopTools:
  T. Hahn, M. Perez-Victoria, Comput. Phys. Commun. 118 (1999) 153
  [hep-ph/9807565]

_EOF_


