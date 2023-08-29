# NexusReader v1.0.0
# Last Modified 8th April 2004
# code (c) Jonathan Jeffery 2002
#
# Reads Nexus format files and puts data into variables
# Then checks for taxonomic equivalence
require 5.0.0;
$VersionTracker = "PerlEQ v1.0";#.12b
use File::Basename;
@SuffixList = ('\..*');
#cheat to get local dir syntax
@LocalPathDetails = fileparse('DoesNotExistAtAll.Hmmm', @SuffixList);
$LocalPath = $LocalPathDetails[1];

#open(TEST,  '>', 'Test.txt');#tracker file for bug hunts

print "\n\n   ?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?\n\n";
print "     PerlEQ v1.0 :  SAFE TAXONOMIC REDUCTION\n\n";
print "   ?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?=?\n\n";

print "       Mark Wilkinson  &  Jonathan Jeffery\n";
print "         NHM, London        Univ. Leiden\n";
print "          (Brains)             (Code)\n\n";
print "            \"Still ain't fancy...\"\n\n";





$Change = ", ";
&NewListSplit;

##################################################
#Startup
# Get input file name and open it
print "\nWhat file do you want to use?\n";
print "Enter the name of a file in the local directory, or a file-path.\n";
print "File must be in Nexus format:\n   ";


while(1){
	$TheFile = <stdin>;
	chomp ($TheFile);
	$TheFile =~ s/\"//g; #get rid of quotes
	
	if((-e $TheFile) == 1){last};
	print "\nCan\'t find that file (wrong spelling or file path?).";
	print "\nPlease try again:\n   ";
	}

@FileDetails = fileparse($TheFile, @SuffixList);
$MyPath = $FileDetails[1];
$InFile = $FileDetails[0] . $FileDetails[2];
open(NEXUS,  '<', $TheFile) or die "The file $TheFile could not be found.\n";
$Assumptfile = $TheFile;

###############################################
#Open output file
&OpenOut;
open(OUTPUT, $OutTypeDat, $PathDataFile);


##############################################################
#Get data from Nexus file
$FileToClean = "NEXUS";
&ImportFile; #output is @OriginalNEXUSstore and @FormatNEXUSstore
@UserNEXUSstore = @OriginalNEXUSstore;
@UserFORMATstore = @FormatNEXUSstore;

#########################################################
#Ask about creating a modified Nexus file
&LayoutCheck;	#returns $Inverted 0/1 and $Interleaf 0/1, $RespectCase;

#Get various parameters of the Nexus file
&NexusCheck; 		#returns $NexusCheck 0/1
&NumberTAXA;		#returns $NTAX, $Newton (newtaxa)
&NumberCHARS;		#returns $NCHAR, $MatrixLabels (left|right|no)
@SymbolsII = ();
&CharSymbols;		#returns @SymbolsII

$Matchchar = "";
$Missing = "";
$Gap = "";
&MissingSymbol;		#returns $Missing, $Gap, $GapMode (missing|newstate) and $Matchchar

@From = ();
@To = ();
&EquateSymbols;		#returns @To and @From equates
&DataTypeCheck;		#returns $DataType (standard, dna, rna, nucleotide, protein)
			# and sets special charsects and equates

#open new nexus file
&NewNexusOpening;

#make sure symbolsII data has no repeated symbols
@tempII = ();
foreach $G (@SymbolsII){
	$X = quotemeta($G);
	if(grep(/$X/, @tempII) == 0){
		push(@tempII, $G);
		}
	}
@SymbolsII = @tempII;

#convert all to upper-case
unless($RespectCase == 1){
	$Missing =~ tr/a-z/A-Z/;
	$Matchchar =~ tr/a-z/A-Z/;
	$Gap =~ tr/a-z/A-Z/;

	@Temp = ();
	@TempA = ();

	for($MO = 0; $MO <= $#To; $MO++){
		$TA = $To[$MO];
		$TB = $From[$MO];
		$TA =~ tr/a-z/A-Z/;
		$TB =~ tr/a-z/A-Z/;
		push(@Temp, $TA);
		push(@TempA, $TB);
		}
	@To = @Temp;
	@From = @TempA;

	@Temp = ();
	foreach $MO (@SymbolsII){
		$MO =~ tr/a-z/A-Z/;
		push(@Temp, $MO);
		}
	@SymbolsII = @Temp;
	}

##################################################
#Create array of ASCII values of Symbols, and
#of 'from' equates for &charcheck
#foreach $i(@SymbolsII){push(@SymbolsListII, ord($i))}
@SymbolsListII = @SymbolsII;
push(@SymbolsListII, $Missing);
if($FoundGap == 1){
	push(@SymbolsListII, $Gap);
	}
#foreach $i(@From){push(@FromII, ord($i))}
##################################################


@TaxLabels = (); # this is every taxon name (or default number)
%TaxLabelTotalArray = (); #this is every taxon name, number and taxset
&GetTaxLabels;
@CharLabels = (); # this is every char name (or default number)
%TotalCharLabels = (); #this is every char name, number and charset
&GetCharLabels;
&GetCharStateLabels;


########################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Read in data Matrix, whether interleaved or not

$AreTherePolym = 0;

#get to matrix
for($StartCol = 0; $StartCol <= $#UserNEXUSstore; $StartCol++){
	$TestLine = $UserNEXUSstore[$StartCol];
	chomp($TestLine);
	if($TestLine =~ /^matrix$/i){
		last;
		}
	}

#set default matrix species labels...
@FromMatrixSpecies = ();
for($bg = 1; $bg <= $NTAX; $bg++){
	push(@FromMatrixSpecies, $bg);
	}
#...and set default matrix char labels...
$FromMatrixCharTest = 0;
@FromMatrixChars = ();
for($bg = 1; $bg <= $NCHAR; $bg++){
	push(@FromMatrixChars, $bg);
	}

if($Inverted == 1){
	$Temp = $NTAX;
	$NTAX = $NCHAR;
	$NCHAR = $Temp;
	}

if($Interleaf == 1){
	@UnrecList = ();
	&ReadInterLeafNexus;
	}

else {
	@UnrecList = ();
	&ReadNexus;
	};

#Convert transposed matrix
if($Inverted == 1){
	$Temp = $NTAX;
	$NTAX = $NCHAR;
	$NCHAR = $Temp;
	if($FromMatrixTest == 1){
		#@FromMatrixChars = @FromMatrixSpecies;
		#$FromMatrixCharTest = 1;
		$chartramp = 0;
		foreach $oq(@FromMatrixChars){
			$chartramp++;
			unless($oq eq "_"){
				$TotalCharLabels{"$oq"} = [('char', $chartramp)]
				}
			}
		}
	@FromMatrixSpecies = @TaxLabels;
	$FromMatrixTest = 0;

	@TempArray = ();

	for($i = 0; $i < $NTAX; $i++){
		@TempTemp = ();
		for($j = 0; $j <= $#NexusArray; $j++){
			$ChTemp = $NexusArray[$j][$i];
			push(@TempTemp, $ChTemp);
			}
		push(@TempArray, [@TempTemp]);
		}
	@NexusArray = @TempArray;
	}

#Put in Taxon labels
#if none on matrix, use tax labels...
if($MatrixLabels eq "no"){
	for($kl = 0; $kl<= $#NexusArray; $kl++){
		@Tempestuous = @{$NexusArray[$kl]};
		$temping = $TaxLabels[$kl];
		unshift(@Tempestuous, $temping);
		$NexusArray[$kl] = [@Tempestuous];
		}
	}

#if labels on matrix...
else{
	for($kl = 0; $kl<= $#NexusArray; $kl++){
		@Tempestuous = @{$NexusArray[$kl]};

		#add non-defaults from taxlabels by preference...
		if($foundTaxLabels == 1){
			$temping = $TaxLabels[$kl];
			}

		#else add non-defaults from matrix...
		elsif(($FromMatrixTest == 1) and ($foundTaxLabels == 0)){
			$temping = $FromMatrixSpecies[$kl];
			$speciestramp = 0;
			foreach $oq(@FromMatrixSpecies){
				$speciestramp++;
				unless($oq eq "_"){
					$TaxLabelTotalArray{"$oq"} = [('taxon', $speciestramp)]
					}
				}
			}

		#else add defaults from taxlabels...
		else{
			$temping = $TaxLabels[$kl];
			}
		unshift(@Tempestuous, $temping);
		$NexusArray[$kl] = [@Tempestuous];
		}
	}

#get max char lengths for printing
@LengthNexusArray = @NexusArray;
&MaxCharLengths;


#Locate Exset and record them
@ExsetNexus = @UserNEXUSstore;
@CharsetOnlyList = (); #3d list of charsets
$AreThereCharSets = 0;
&FindCharacterSets;
&PAUPBLOCK; #read in paup block for later use by &ExcludeFinder and &ExcludeTaxa

# Get big list of all assumptions data
#put here to alow codon sets to be used in exsets
$Assumptfile = $TheFile;
close(NEXUS);
@Typesets = ();
@CodonSets = ();
&AssumptSub; 			#looks for charsets, codonpossets and ctypes
if($AreThereCodons == 1){	#finds and reads codonpossets.
	&CodonReader;
	}
#make default Exset (all included)
for($IU = 0; $IU < $NCHAR; $IU++){push(@FinalExSet, 0)}
&ExcludeFinder;

#Locate Taxsets and  record them
@TaxSetFinalList = (); #3d list of all taxsets

if($SomeExcluded == 0){
	&ExSetFinder;
	}
@FinalNameList = ();
$AreThereTaxSets = 0;
&FindTaxonSets;

#Locate delete statements and record them
#plus how to treat polymorphic taxa (also in PAUP block)
$MStaxa = "variable";
$FoundMStaxa = 0;
&ExcludeTaxa;

#########################################
# Print data collected on size and type of file
$PathInFile = $InFile;
unless($MyPath eq $LocalPath){
	$PathInFile = $MyPath . $InFile;
	}

#Make header
print OUTPUT "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\"
\"http://www.w3.org/TR/html4/loose.dtd\">
<html lang=\"en\">
<head>
<title>Data from $PathInFile checked for equivalence by $VersionTracker</title>
<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\">
<style>
body, p	{font-family:\"Times New Roman\", Times, serif;
         font-weight:normal;
         font-style:normal;
         font-size:12pt;
         color:#000000;
         text-decoration:none;}

body {background-color:#FFFFFF;}

p, p.Grey, p.GreyNex {display:inline;}
p.Grey, p.GreyNex {color:#CCCCCC;}
p.GreyNex {font-family:\"Courier New\", Courier, mono;}

h1,h2,h3 {font-weight:bold;
          font-size:18pt;
          }
h2 {font-size:16pt;
    text-decoration:underline;}
h3 {font-size:12pt;
    text-decoration:underline;}

table, table.top {width:100%; border-collapse:collapse; white-space:nowrap}
table.top {background:#52C7ED;}
table.nexus {width:0;}

tr, tr.clear {border-spacing:0px;
              border-width:0px;
              padding:0px;
              white-space:nowrap;}
tr.clear {background:#F7F7F7;}

th, th.top {font-weight:bold;
            font-size:24pt;
            padding:5px;}
th.Next	{font-size:16pt;}

td, td.nexusSp, td.nexusDat {
		border-width:0px;
		padding:0px;
		white-space:nowrap;}
td.nexusEq {padding:5px;}
td.nexusDat {font-family : \"Courier New\", Courier, mono;}
</style>
</head>\n\n
<body>";

#title table
print OUTPUT "\n<table class=\"top\">\n\t<tr>";
print OUTPUT "\n\t<th class=\"top\"><a name=\"Top\">Data from <i>$PathInFile</i></a></th>
\t</tr>
\t<tr>
\t\t<th class=\"Next\">Checked for equivalence by <i>$VersionTracker</i>";
&CreateDate;
print OUTPUT "\n\t\t<br>Created $Days[$CreateDay], $CreateDate $Months[$CreateMonth] $CreateYear at ";
printf OUTPUT ("%02d:%02d:%02d</th></tr></table>", $CreateHour, $CreateMin, $CreateSec);
print OUTPUT "\n<p><br>Short-cut to final <A HREF=\"$DataFile\#STRList\">deletion-list</A>\n</p>";


#check ndel & nexclud
$Dotter = 0;
if($SomeDeleted == 1){
	for($DeleteThem = 0; $DeleteThem <= $#NexusArray; $DeleteThem++){
		if($DeletionFinalList[$DeleteThem] == 0){
			$Dotter++;
			}
		}
	}
$Totter = 0;
for($DD = 0; $DD <= $#FinalExSet; $DD++){
	if($FinalExSet[$DD] == 0){next}
	else {$Totter++}
	}


print OUTPUT "\n<hr>\n";
print OUTPUT  "<h2>1. Details of Input File</h2>";
if($Dotter > 0){
	$DelBoy = "\(of which $Dotter deleted\)";
	}
print OUTPUT "<p>Number of taxa = <b>$NTAX $DelBoy</b>\n";
if($Totter > 0){
	$DelBoy = "\(of which $Totter excluded\)";
	}
print OUTPUT "<br>Number of characters = <b>$NCHAR $DelBoy</b>\n";
print OUTPUT "<br>Symbols = <b>@SymbolsII</b>\n";
print OUTPUT "<br>Missing = <b>$Missing</b>\n";
if($FoundGap == 1){
	print OUTPUT "<br>Gap = <b>$Gap</b> \(treated as $GapMode\)\n";
	}

if($AreTherePolym == 1){
	print OUTPUT "<br>Multistate scores treated as $MStaxa";
	if($MStaxa eq "variable"){
		print OUTPUT ": \(...\)=polymorphic, \{...\}=uncertain";
		}
	if($PolyKill == 1){
		print OUTPUT "<br><i>Note that polymorphic scores are, for the moment, recoded as \'$Missing\'<\/i>";
		}
	print OUTPUT "<br>\n";


	}
print OUTPUT "\n<br></p>\n\n";

#list of excluded chars and taxa
unless(@From == 0){
	print OUTPUT "<p><br>\n";
	for($i = 0; $i <= $#From; $i++){
		print OUTPUT "Character code \'$From[$i]\' changed to \'$To[$i]\'<br>\n";
		}
	print OUTPUT "<br></p>\n";
	}

$ExcludeCount = 0;
for($DD = 0; $DD <= $#FinalExSet; $DD++){
	if($FinalExSet[$DD] == 0){next}
	else {
		$ExcludeCount = 1;
		}
	}

$EXCLUDEDcountTotal = 0;
if($ExcludeCount != 0){
	$Plup = "";
	if($Totter > 1){$Plup = "s"}
	if($SomeExcluded == 1){
		print OUTPUT "\n<ul><li>$Totter character$Plup excluded via a <i>PAUP</i>-block statement \(printed in pale grey\):</li>\n";
		$Plup = "";
		}
	else{
		$Quotey = "";
		unless($UsedExset =~ /^\'|\'$/){
			$Quotey = "\'";
			}
		if(length($UsedExset) == 0){
			print OUTPUT "<ul><li>$Totter character$Plup excluded via an \'eliminate\' statement \(printed in pale grey\):</li>\n";
			$Plup = "";
			}
		else{print OUTPUT "<ul><li>$Quotey$UsedExset$Quotey used as character exclusion set \(printed in pale grey\):</li>\n";
			$Plup = " excluded";
			}
		}
	$AddNote = 0;

	print OUTPUT "\n<blockquote>";
	$CommaRight = "";
	$printstring = "";
	$printstleng = 0;

	for($DD = 0; $DD <= $#FinalExSet; $DD++){
		if($FinalExSet[$DD] == 0){next}
		else {
			$EXCLUDEDcountTotal++;
			$Asterix = "";
			$ExcludeCount = 1;
			foreach $i (@EliminateDiff){
				if(($i == ($DD + 1)) and (length($UsedExset) != 0)){$Asterix = "\*";
					$AddNote = 1;
					}
				}

			if($printstleng >= 60){
				$CommaRight = ",<br>";
				$printstleng =  0;
				}
			$printstring .= "$CommaRight". ($DD + 1) . "$Asterix";
			$printstleng = $printstleng + length("$CommaRight". ($DD + 1) . "$Asterix");
			$CommaRight = ", ";

			}
		}
	print OUTPUT "\n$printstring$Plup</blockquote>";
	if($AddNote == 1){
		print OUTPUT "<p><b>*</b> = Character excluded by an \'eliminate\' statement<br></p>\n"
		}
	}

elsif($ExcludeCount == 0){
	print OUTPUT "<ul><li>All characters included in analysis</li>\n\n";
	}

if(@UnrecList > 0){
	for($SDRR = 0; $SDRR <= $#UnrecList; $SDRR++){
		print OUTPUT "<li>$UnrecList[$SDRR]</li>\n";
		}
	}

#delete taxa now and record it
$DELETEDcount = 0;
if($SomeDeleted == 1){
	$Plup = "on";
	if($Dotter > 1){$Plup = "a"}
	print OUTPUT "\n<li>$Dotter tax$Plup deleted by user \(printed in pale grey\):</li>\n<blockquote>";
	$CommaRight = "";
	$printstring = "";
	$printstleng = 0;

	for($DeleteThem = 0; $DeleteThem <= $#NexusArray; $DeleteThem++){
		if($DeletionFinalList[$DeleteThem] == 0){
			$DELETEDcount++;
			$quoteless = $NexusArray[$DeleteThem][0];
			$quoteless =~ s/^\'//;
			$quoteless =~ s/\'$//;

			if($printstleng >= 60){
				$CommaRight = ",<br>";
				$printstleng =  0;
				}
			$printstring .= "$CommaRight<i>$quoteless<\/i>";
			$printstleng = $printstleng + length("$CommaRight" . "$quoteless");
			$CommaRight = ", ";

			}
		}
	print OUTPUT "$printstring</blockquote><br>\n";
	}

print OUTPUT "</ul>";

#Print Nexus matrix
print OUTPUT "\n<hr>\n";
 print OUTPUT  "<h2>2. Original Data Matrix</h2>";


$CharTaxCol = "Taxa ";
$CharTaxRow = "Chars";
@NexusArrayINTER = @NexusArray;



$PrintTo = "OUTPUT";
$Explain = "<p><br><b>The matrix is large; for clarity it is printed in blocks of 50 Characters</b></p>\n";
$BrakishO = " ";
$BrakishC = " ";
&INTERLEAFPRINT;

print OUTPUT "\n\n<br><hr>\n";
print OUTPUT  "<h2>3. Data on Individual Characters</h2>";


# looking for the assumpts to be used
# Get big list of all assumptions data
&DefAssumptSub; 	#gets default chartypes
$Misunderstood = 0;
&SettingType;		#sets default type

if($CTypeCheck > 0){&CtypeReadin}
if($CTypeChanges == 0){&ExaminationTypeset}

@StartTime = localtime(time);
@T1 = localtime(time);

print "\n\n???????????????????????????????????????????????????????????????????????????\n";
print   "\n Running analysis (this may take quite a while for really large data sets)\n";
print   "\n???????????????????????????????????????????????????????????????????????????\n";
print   "\nPercent Completed:\n";

################################################

#Change multistate ordered to a range of value
#(I think this is the incorrect way to handle such characters
#but PAUP* does it this way...)
@RecodeArray = @NexusArray;
@RecodeSymbols = @SymbolsII;
if(($FoundGap == 1) and ($GapMode ne "missing")){
	push(@RecodeSymbols, $Gap);
	}

&OrderRecode;
@NexusArray = @RecodeArray;

################################################
#Get data on Informativeness
@NexusArraySub = @NexusArray;
&Informative;

@CharStates = @CharStatesSub;
@StateList = @StateListSub; # 3D list, contains numbers of states for each char
@PolyList = @PolyListSub; # 3D list contains all polymorphic codes for each char

################################################
#print a list of char types

print OUTPUT "\n<table class=\"nexus\">
<tr>
\t<td class=\"nexusSp\">&nbsp;</td>
\t<td class=\"nexusSp\">&nbsp;</td>
\t<td class=\"nexusSp\"><div align=\"left\"><b>States:</b></div></td>
</tr>";

print OUTPUT "\n<td class=\"nexusSp\"><div align=\"center\"><b>Character&nbsp;&nbsp;</b></div></td>\n\t<td class=\"nexusSp\"><b>Type</b></td> ";

for($v = 0; $v <= $#SymbolsII; $v++){
	print OUTPUT "\n\t<td class=\"nexusSp\"><div align=\"right\"><b>&nbsp;&nbsp;&nbsp;$SymbolsII[$v]</b></div></td>"
	}
if($FoundGap == 1){
	print OUTPUT "\n\t<td class=\"nexusSp\"><b>&nbsp;&nbsp;Gap&nbsp;&nbsp;</b></td>";
	}

if($FoundPoly == 1){
	print OUTPUT "\n\t<td class=\"nexusSp\"><b>&nbsp;&nbsp;MultiSt&nbsp;&nbsp;</b></td>";
	}

if(($FoundGap == 0) and ($FoundPoly == 0)){print OUTPUT "\n\t<td class=\"nexusSp\"><b>&nbsp;&nbsp;Missing&nbsp;&nbsp;</b></td>\n</tr>";}
else{print OUTPUT "\n\t<td class=\"nexusSp\"><b>Missing</b></td>\n</tr>";}

$AllUnderstood = 0;
for($u = 0; $u < $NCHAR; $u++){
	print OUTPUT "\n<tr>\n\t<td class=\"nexusSp\"><div align=\"center\">";
	if($FinalExSet[$u] == 1){
		print OUTPUT "<p class=\"Grey\">";
		}
	print OUTPUT $u + 1;
	if($FinalExSet[$u] == 1){
			print OUTPUT "</p>";
		}
	print OUTPUT "</div></td>";
	if($FinalExSet[$u] == 1){
		print OUTPUT "\n\t<td class=\"nexusSp\"><p class=\"Grey\"><i>excluded&nbsp;</i><\/p></td>\n</tr>";
		next;
		}
	else{
		print OUTPUT "\n\t<td class=\"nexusSp\">$CharTypes2[$u]\&nbsp;";


		if($Understood[$u] == 0){
			print OUTPUT "</td>"
			}
		else{
			$temp2--;
			print OUTPUT "\*</td>";
			$AllUnderstood = 1
			}
		@temp = @{$CharStates[$u]};
		@poly = @{$PolyList[$u]};
		$NoPoly = @poly;

		#print charstates and gap
		for($v = 0; $v < $#temp; $v++){

			print OUTPUT "\n\t<td class=\"nexusSp\">";
			if(($FoundGap == 1) and ($v == ($#temp - 1))){
				print OUTPUT "<div align=\"center\">$temp[$v]";
				}
			else{print OUTPUT "<div align=\"right\">&nbsp;&nbsp;&nbsp;$temp[$v]";}
			print OUTPUT "</div></td>";
			}

		#print polymorphic
		if($FoundPoly == 1){
			print OUTPUT "\n\t<td class=\"nexusSp\"><div align=\"center\">$NoPoly</div></td>";
			}

		#print missing
		print OUTPUT "\n\t<td class=\"nexusSp\"><div align=\"center\">" . ($temp[$#temp] - $Dotter) . "</div></td>\n</tr>";
		}
	}
if($AllUnderstood == 0){print OUTPUT "</table><p><br><br></p>"}
else{print OUTPUT "</Table><p><br><b>*</b> = Original type not understood, character set to unordered<br><br></p>\n\n"}


#############################
#Print data on distribution of missing data
$WhatCounting = "";
if(($DELETEDcount > 0) or ($EXCLUDEDcountTotal > 0)){
	if(($DELETEDcount > 0) and ($EXCLUDEDcountTotal > 0)){
		$WhatCounting = "\n \(not counting deleted taxa or excluded characters\)";
		}
	if(($DELETEDcount > 0) and ($EXCLUDEDcountTotal == 0)){
		$WhatCounting = " \(not counting excluded characters\)";
		}
	if(($DELETEDcount == 0) and ($EXCLUDEDcountTotal > 0)){
		$WhatCounting = " \(not counting deleted taxa\)";
		}
	}

print OUTPUT "\n<h3>3a. Statistics on Missing Data$WhatCounting</h3>\n";

for($u = 0; $u <= ($NTAX -1); $u++){
	push(@DistribMissing, 0)
	}
$TotalMissing = 0;

for($u = 0; $u < $NCHAR; $u++){
	if($FinalExSet[$u] == 1){
		next
		}
	else{
		@temp = @{$CharStates[$u]};
		if(($FoundGap == 1) and ($GapMode eq "missing")){
			$y = $temp[$#temp] + $temp[$#temp - 1] - $DELETEDcount;
			}
		else{$y = $temp[$#temp] - $DELETEDcount}

		$TotalMissing = $TotalMissing + $y;
		$DistribMissing[$y]++;
		}
	}
print OUTPUT "<table class=\"nexus\">\n";
for($u = 0; $u <= ($#DistribMissing); $u++){
	if($u == 1){$Entry = "try"}
	else{$Entry = "tries"}
	if($DistribMissing[$u] == 0){next}
	elsif($DistribMissing[$u] == 1){$Character = "ter has"}
	else{$Character = "ters have"}
	print OUTPUT "\t<tr align=\"center\">\n\t\t<td class=\"nexusSp\"><div align=\"right\"><b>$DistribMissing[$u]&nbsp;</b></div></td>\n\t\t<td class=\"nexusSp\">charac$Character got</td>\n\t\t<td class=\"nexusSp\"><b>&nbsp;$u&nbsp;</b></td>\n\t\t<td class=\"nexusSp\">missing en$Entry</td>\n</tr>\n";
	}
print OUTPUT "</table>\n";

$TotalChars = ($NTAX - $DELETEDcount) * ($NCHAR - $EXCLUDEDcountTotal);
$Proportion = ($TotalMissing/$TotalChars)  * 100;
$MeanPerChar = $TotalMissing/($NCHAR - $EXCLUDEDcountTotal);
$MeanPerTax = $TotalMissing/($NTAX - $DELETEDcount);
print OUTPUT "\n<p><br>Total number of entries = <b>$TotalChars</b>";
if($TotalMissing > 0){
	print OUTPUT "\n<br>Total number of missing entries = <b>$TotalMissing \(";
	printf OUTPUT ("%.2f",$Proportion);
	print OUTPUT "\%\)</b>";
	print OUTPUT "\n<Table class=\"nexus\">\n<tr>\n\t<td class=\"nexusSp\">Mean number of missing entries =</td>\n\t<td  class=\"nexusSp\"><b>";
	printf OUTPUT (" %.2f", $MeanPerTax);
	print OUTPUT " per taxon</b></td>\n\t</tr>\n<tr>\n\t<td class=\"nexusSp\">&nbsp;</td>\n\t<td class=\"nexusSp\">";
	printf OUTPUT ("<b> %.2f", $MeanPerChar);
	print OUTPUT " per character</b></td>\n</tr>\n</table></p>";
	}
else{print OUTPUT "\n<br><b>No missing entries</b></p>";}
print OUTPUT "\n<br>";

@MissingNexusArray = @NexusArray;
$TotalInvalid = 0;
$Excess = 0;
$TimeIt = 0;
$OldTotUp = 0;
&Comparisons; #outputs 3D array @FinalEQdata [name, spacingNo, Missing, Equivs]
&SuffixFixer;


if($TotalInvalid >0){$PropInval = ($TotalInvalid/$TestTotalComp) * 100}
$TotalInvalidStore = $TotalInvalid;

print OUTPUT "\n<h3>3b. Statistics on Invalid Comparisons due to Missing Data$WhatCounting</h3>";
print OUTPUT "\n<p>Total number of comparisons = <b>$TestTotalComp</b>";
if($TotalInvalid > 0){
	print OUTPUT "\n<br>Total number of invalid comparisons = <b>$TotalInvalid \(";
	printf OUTPUT ("%.2f",$PropInval);
	print OUTPUT "\%\)</b></p>";
	}
else{print OUTPUT  "\n<br><b>No invalid comparisons due to missing entries</b></p>"}

###############################################
#Compare missing data

print OUTPUT "\n<p><br><br><hr></p>";
print OUTPUT "\n<h2>4. Equivalents Considering Missing Entries Only</h2>";
print OUTPUT "\n\n<Table class=\"nexus\">\n<tr>\n\t<td class=\"nexusEq\"><b>Index Taxon</b></td>";
print OUTPUT "\n\t<td class=\"nexusEq\"><div align=\"center\"><b>\[Missing\/Uncertain\]</b></div></td>\n\t<td class=\"nexusEq\"><b>Equivalents</b></td>\n</tr>";
&PrintEQ;
print OUTPUT "\n<br></p>";

###############################################
#Change dataset to get rid of uninformative chars  @CharTypes2 @CharStates = states
@ChangedData = ();
$TestChangedData = 0;
$HasAnyThingChanged = 0;

@UninfNexusArray = @NexusArray;
@CharStatesInf = @CharStates;
print OUTPUT "\n";

$DoAgain = 0;

for($i = 0; $i < $NCHAR; $i++){
	$DoAgain = 0;


	#get data on the states of this char
	$t = $i + 1;
	$states = $#SymbolsII + 1;
	@CharStatesSub = ();
	@StateListSub = ();

	for($u = 0; $u <= $states; $u++){
		push(@StateListSub, 0)
		}
	if($FoundGap == 1){
		push(@StateListSub, '0');
		}
	&LittleInf; #generates @StateListSub and @TempPoly

	@TempStateData = @StateListSub;
	@TempPolyData = @TempPoly; # these used in &OverallInform & &Redistrib

	#go to appropriate testing-sub
	if($CharTypes2[$i] =~ /^unord/){
		&ReplaceUnord
		}
	elsif(($CharTypes2[$i] =~ /^ord/) or ($CharTypes2[$i] =~ /dollo\.up/)or ($CharTypes2[$i] =~ /dollo\.dn/)){
		&ReplaceLower;
		&ReplaceUpper;
		}
	elsif($CharTypes2[$i] =~ /irrev\.up/){&ReplaceUpper}
	elsif($CharTypes2[$i] =~ /irrev\.dn/){&ReplaceLower}

	#if changes have been made, try char again
	if($DoAgain == 1){
		$DoAgain = 0;
		$i--;
		}
	}

#take nexus stripped of uniformative chars, and
#puts it into the input nexus for the informative sub
@NexusArraySub = @UninfNexusArray;
&Informative;

#from informative sub, a 3D list of each character,
#and the number of entries it has at each state
@CharStatesInf = @CharStatesSub;
$BadYes = 0;
@Badlist = ();

@InfSymbols2 = @SymbolsII;
if(($FoundGap == 1) and ($GapMode ne "missing")){
	push(@InfSymbols2, $Gap);
	}

&OverallInform;

if($TestChangedData == 0){print OUTPUT "\n\n<p><b>No uniformative data in Nexus</b><br></p>"}

else{
	print OUTPUT "\n<hr><h2>5. Recoding Uninformative Data</h2>";

	#List Changed data;
	if($HasAnyThingChanged == 1){
		print OUTPUT "\n<h3>5a. Uninformative Codings</h3>";
		print OUTPUT "\n<table class=\"nexus\">\n<tr>\n\t<td class=\"nexusSp\"><b>Taxon&nbsp;</b></td> ";
		print OUTPUT "\n\t<td class=\"nexusSp\"><b>&nbsp;Character&nbsp;</b></td>\n\t<td class=\"nexusSp\"><b>&nbsp;State</b></td>\n</tr>";
	#loop to get changes printed in order
		@TempSymbols2 = @SymbolsII;
		if(($FoundGap == 1) and ($GapMode ne "missing")){
			push(@TempSymbols2, $Gap);
			}

		for($h = 0; $h <= $NTAX; $h++){
			for($w = 0; $w <= $NCHAR; $w++){
				for($f = 0; $f <= $#TempSymbols2; $f++){
					for($g = 0; $g <= $#ChangedData; $g++){
						if(
							($ChangedData[$g][0] == $h) and
							($ChangedData[$g][1] == $w) and
							($ChangedData[$g][2] eq $TempSymbols2[$f])
							){
							&PrintChange;
							}
						else{next}
						}
					}
				}
			}
		print OUTPUT "\n</table>\n"
		}#end of if

	#print out chars reduced to all-missing

	if($BadYes == 1){
		print OUTPUT "\n<h3>5b. Characters With No Parsimony-Informative Data (replaced with \'$Missing\')</h3>\n<p><blockquote>";
		$printout = "";
		$Stringlength = 0;

		#make list with no excluded chars...
		@PrintBad = ();
		for($TT = 0; $TT <= $#Badlist; $TT++){
			$CharToCheck = $Badlist[$TT];
			if($FinalExSet[$CharToCheck - 1] == 0){
				push(@PrintBad, $CharToCheck);
				}
			}


		#prints out a compressed list, with 60chr line-length
		for($TT = 0; $TT <= $#PrintBad ; $TT++){
			$printadd = "";
			$This = $PrintBad [$TT];
			$That = $PrintBad [$TT + 2];
			$Adder = 2;

			if($That == ($This + $Adder)){
				$printadd =  " $This-";
				$Adder++;

				$TheOther = $PrintBad[$TT + $Adder];
				until($TheOther !=  ($This + $Adder)){
					$That = $TheOther;
					$Adder++;
					$TheOther = $PrintBad[$TT + $Adder];
					}
				$printadd .= "$That,";
				$TT = $TT + ($Adder - 1);
				}
			else{
				$printadd = " $This,";
				}
			$Stringlength = $Stringlength + length($printadd);

			if($Stringlength > 70){
				$printout .= "<br>\n\t";
				$printadd =~ s/^ //;
				$Stringlength = length($printadd);

				}
			$printout .= $printadd;
			}

		chop($printout);
		$printout =~ s/^ //;
		print OUTPUT "$printout</blockquote></p>\n";
		}
	else{print OUTPUT "\n<p><br></p>"}

	@NexusArraySub = @UninfNexusArray;
	&Informative;

	#list of no of states of each character, after uninf. stripped
	@CharStatesUninf = @CharStatesSub;
	@StateList = @StateListSub;

	#############################
	#Print data on distribution of uninformative data
	print OUTPUT "\n<hr>\n<h2>6. Statistics on Uninformative Data$WhatCounting</h2>\n";

	for($u = 0; $u < $NTAX; $u++){
		push(@DistribUninf, 0)
		}
	$TotalUninf = 0;
	$GrandTotalMU = 0;
	for($u = 0; $u < $NCHAR; $u++){
		if($FinalExSet[$u] == 1){
			next
			}
		else{
			#get data for informative chars
			@tempUnif = @{$CharStatesUninf[$u]};
			@tempMiss = @{$CharStates[$u]};
			if(($FoundGap == 1) and ($GapMode eq "missing")){
				$Difference = $tempUnif[$#tempUnif] - $tempMiss[$#tempUnif] + $tempUnif[$#tempUnif - 1] - $tempMiss[$#tempUnif - 1];
				$GrandTotalMU = $GrandTotalMU + $tempUnif[$#tempUnif] + $tempUnif[$#tempUnif - 1]  - $DELETEDcount;
				}
			else{
				$Difference = $tempUnif[$#tempUnif] - $tempMiss[$#tempUnif];
				$GrandTotalMU = $GrandTotalMU + $tempUnif[$#tempUnif] - $DELETEDcount;
				}
			$TotalUninf = $TotalUninf + $Difference;
			$DistribUninf[$Difference]++;
			}
		}
	print OUTPUT "<table class=\"nexus\">\n";
	for($u = 0; $u <= $#DistribUninf; $u++){
		if($u == 1){$Entry = "try"}
		else{$Entry = "tries"}
		if($DistribUninf[$u] == 0){next}
		elsif($DistribUninf[$u] == 1){$Character = "ter has"}
		else{$Character = "ters have"}
		print OUTPUT "\n<tr align=\"center\">\n\t<td class=\"nexusSp\"><div align=\"right\"><b>$DistribUninf[$u]&nbsp;</b></div></td>";
		print OUTPUT "\n\t<td class=\"nexusSp\">charac$Character got</td>\n\t<td class=\"nexusSp\"><b>&nbsp;$u&nbsp;</b></td>\n\t<td class=\"nexusSp\">scored, but uninformative, en$Entry</td>\n</tr>\n";
		}
	print OUTPUT "</table>\n";

	$TotalChars = ($NTAX - $DELETEDcount) * ($NCHAR - $EXCLUDEDcountTotal);
	$Proportion = ($TotalUninf/$TotalChars) * 100;
	$MeanPerChar = $TotalUninf/($NCHAR - $EXCLUDEDcountTotal);
	$MeanPerTax = $TotalUninf/($NTAX - $DELETEDcount);
	print OUTPUT "\n<p><br>Total number of entries =<b> $TotalChars</b>";
	print OUTPUT "\n<br>Total number of uninformative entries = <b>$TotalUninf \(";
	printf OUTPUT ("%.2f",$Proportion);
	print OUTPUT "\%\)</b></p>";
	print OUTPUT "\n<table class=\"nexus\">\n<tr>\n\t<td class=\"nexusSp\">Mean number of uninformative entries =</td>\n\t<td class=\"nexusSp\">";
	printf OUTPUT ("<b> %.2f", $MeanPerTax);
	print OUTPUT " per taxon</b></td>\n</tr>\n<tr>\n\t<td class=\"nexusSp\">&nbsp;</td>\n\t<td class=\"nexusSp\">";
	printf OUTPUT ("<b> %.2f", $MeanPerChar);
	print OUTPUT " per character</b></td>\n</tr>\n</table>\n";

	@MissingNexusArray = @UninfNexusArray;
	$TotalInvalid = 0;
	&Comparisons;
	&SuffixFixer;
	$UninfInval = $TotalInvalid - $TotalInvalidStore;
	$PropInval = $UninfInval/$TestTotalComp;

	print OUTPUT "<h3><br>6a. Statistics on Invalid Comparisons due to Uninformative Data$WhatCounting</h3>\n";
	print OUTPUT "\n<p>Total number of comparisons = <b>$TestTotalComp</b>";
	print OUTPUT "\n<br>Total number of invalid comparisons = <b>$UninfInval \(";
	printf OUTPUT ("%.2f",$PropInval);
	print OUTPUT "\%\)</b><br></p>";
	print OUTPUT "<h3><br>6b. Combined statistics on Missing <i>and</i> Uninformative Data$WhatCounting</h3>\n";
	$PropInval = ($TotalInvalid/$TestTotalComp) * 100;
	$Proportion = ($GrandTotalMU/$TotalChars) * 100;
	$MeanPerChar = $GrandTotalMU/($NCHAR - $EXCLUDEDcountTotal);
	$MeanPerTax = $GrandTotalMU/($NTAX - $DELETEDcount);

	print OUTPUT "\n<p>Total number of entries = <b>$TotalChars</b>";
	print OUTPUT "\n<br>Total number of missing and uninformative entries = <b>$GrandTotalMU (";
	printf OUTPUT ("%.2f",$Proportion);
	print OUTPUT "\%\)</b></p>";
	print OUTPUT "\n<table class=\"nexus\">\n<tr>\n\t<td class=\"nexusSp\">Mean number of missing and uninformative entries =</td>\n\t<td class=\"nexusSp\">";
	printf OUTPUT ("<b> %.2f", $MeanPerTax);
	print OUTPUT " per taxon</b></td>\n</tr>\n<tr>\n\t<td class=\"nexusSp\">&nbsp;</td>\n\t<td class=\"nexusSp\">";
	printf OUTPUT ("<b> %.2f", $MeanPerChar);
	print OUTPUT " per character</b></td>\n</tr>\n</table>\n";
	print OUTPUT "\n<p>Total number of comparisons = <b>$TestTotalComp</b>";
	print OUTPUT "\n<br>Total number of invalid comparisons = <b>$TotalInvalid \(";
	printf OUTPUT ("%.2f",$PropInval);
	print OUTPUT "\%\)</b><br><br></p>";

	##############################
	#Print data on Tax Equiv after unif chars have been removed
	print OUTPUT "<hr><h2>7. Equivalents Considering Missing <i>and</i> Uninformative Entries</h2>";
	print OUTPUT "\n\n<table class=\"nexus\">\n<tr>\n\t<td class=\"nexusEq\"><b>Index Taxon</b></td>";
	print OUTPUT "\n\t<td class=\"nexusEq\"><div align = \"center\"><b>\[Missing\/Uncertain\]</b></div></td>\n\t<td class=\"nexusEq\"><b>Equivalents</b></td>\n</tr>";
	&PrintEQ;
	}#end of else to print uninf data

#print key to symbols
print OUTPUT "\n</p>\n<hr>";
print OUTPUT "\n<h2>Key to Taxonomic Equivalence</h2>";
print OUTPUT "\n<table class=\"nexus\">\n<tr>\n\t<td class=\"nexusSp\" colspan=\"2\"><b>Equivalence Scores:</b></td>\n</tr>";
print OUTPUT "\n<tr class=\"clear\">\n\t<td class=\"nexusSp\"><b>A</b> - </td>\n\t<td class=\"nexusSp\">Actual equivalents (symmetric)</td>\n</tr>";
print OUTPUT "\n<tr class=\"clear\">\n\t<td class=\"nexusSp\">&nbsp;</td>\n\t<td class=\"nexusSp\">Equivalent can be safely deleted if Index is retained</td>\n</tr>";
print OUTPUT "\n<tr>\n\t<td class=\"nexusSp\"><b>B</b> - </td>\n\t<td class=\"nexusSp\">Potential equivalents (symmetric)</td>\n</tr>";
print OUTPUT "\n<tr>\n\t<td class=\"nexusSp\">&nbsp;</td>\n\t<td class=\"nexusSp\">Equivalent can be safely deleted if Index is retained</td>\n</tr>";
print OUTPUT "\n<tr class=\"clear\">\n\t<td class=\"nexusSp\"><b>C</b> - </td>\n\t<td class=\"nexusSp\">Potential equivalents (asymmetric all one way)</td>\n</tr>";
print OUTPUT "\n<tr class=\"clear\">\n\t<td class=\"nexusSp\">&nbsp;</td>\n\t<td class=\"nexusSp\">Equivalent can be safely deleted if Index is retained</td>\n</tr>";
print OUTPUT "\n<tr>\n\t<td class=\"nexusSp\"><b>D</b> - </td>\n\t<td class=\"nexusSp\">Potential equivalents (asymmetric both ways)</td>\n</tr>";
print OUTPUT "\n<tr>\n\t<td class=\"nexusSp\">&nbsp;</td>\n\t<td class=\"nexusSp\">Cannot be safely deleted</td>\n</tr>";
print OUTPUT "\n<tr class=\"clear\">\n\t<td class=\"nexusSp\"><b>E</b> - </td>\n\t<td class=\"nexusSp\">Potential equivalents (asymmetric all one way)</td>\n</tr>";
print OUTPUT "\n<tr class=\"clear\">\n\t<td class=\"nexusSp\">&nbsp;</td>\n\t<td class=\"nexusSp\">Index can be safely deleted if equivalent is retained</td>\n</tr>";
print OUTPUT "\n<tr>\n\t<td class=\"nexusSp\" colspan=\"2\"><br><b>Suffixes:</b></td>\n</tr>";
print OUTPUT "\n<tr class=\"clear\">\n\t<td class=\"nexusSp\"><b>*</b> - </td>\n\t<td class=\"nexusSp\">Equivalent taxon must originate from same node as Index taxon in</td>\n</tr>";
print OUTPUT "\n<tr class=\"clear\">\n\t<td class=\"nexusSp\">&nbsp;</td>\n\t<td class=\"nexusSp\">any MPT (assuming no arbitrary resolutions)</td>\n</tr>";
print OUTPUT "\n<tr>\n\t<td class=\"nexusSp\"><b>!</b> - </td>\n\t<td class=\"nexusSp\">Equivalent taxon has no informative data</td></tr>";
print OUTPUT "\n</table>\n\n";


#print suggested deletioin list
print OUTPUT "<p><br></p><hr>\n";
print OUTPUT "<h2><a name=\"STRList\">8. Scope for Safe Taxonomic Reduction</h2></a>\n";
print OUTPUT "<p><b>Return to <A HREF=\"$DataFile\#Top\">top</a></b><br><br>";

@Tit = sort{&TitSorter}(keys(%SuggestDel));


$Is_there_anything = 0;
$Then_how_many = 0;

for($RTV = 0; $RTV <= $#Tit; $RTV++){
	$NaNa = $SuggestDel{"$Tit[$RTV]"};
	if($NaNa == 0){
		$Is_there_anything = 1;
		$Then_how_many++;
		}
	}
if($Is_there_anything == 0){
	print OUTPUT "Sorry, no scope for safe taxonomic reduction\!</p>\n\n";
	}
else{
	$AreThereTaxSets = 1;
	$CorrectPleural = "taxon is";
	if($Then_how_many > 1){$CorrectPleural = "taxa are";}
	print OUTPUT "We suggest the following $CorrectPleural deleted:<blockquote>";
	for($RTV = 0; $RTV <= $#Tit; $RTV++){
		$NaNa = $SuggestDel{"$Tit[$RTV]"};
		if($NaNa == 0){
			$quoteless = $Tit[$RTV];
			$quoteless =~ s/^\'//;
			$quoteless =~ s/\'$//;
			print OUTPUT "\n<i>$quoteless</i><br>";
			}
		}
	print OUTPUT "\n</blockquote>";
	if(($SomeDeleted == 1) and ($Is_there_anything > 0)){
		print OUTPUT "\nThis is in addition to those already deleted by the user<br>";
		}

	print OUTPUT "\nA \'delete\' command has already been added to \'<a href=\"$NexusFile\" target=\"_blank\">$NexusFile</a>\'";
	print OUTPUT "\n\n";
	}


#print new nexus

$CharTaxCol = "Taxa ";
$CharTaxRow = "Chars";

$Change = " ";
&NewListSplit;
print NEWNEXUS "#NEXUS\n";
print NEWNEXUS "\[\!Data from $PathInFile checked for Taxonomic Equivalence";
print NEWNEXUS "\nCreated $Days[$CreateDay], $CreateDate $Months[$CreateMonth] $CreateYear at ";
printf NEWNEXUS ("%02d:%02d:%02d", $CreateHour, $CreateMin, $CreateSec);
print NEWNEXUS "\]";
#get taxon names to use
@NexusTaxonSet = ();
for($khg = 0; $khg <= $#UninfNexusArray; $khg++){
	push(@NexusTaxonSet, $UninfNexusArray[$khg][0]);
	}


#print taxa set
print NEWNEXUS "\n\nBegin Taxa\;\n";
print NEWNEXUS "\tDimensions ntax=$NTAX\;";
print NEWNEXUS "\n\tTaxlabels\n\t";

$LineTracker = 1;
$printout = "";
$Stringlength = 0;

for($khg = 1; ($khg - 1) <= $#NexusTaxonSet; $khg++){
	$Stringlength = $Stringlength + 1 + length($NexusTaxonSet[$khg - 1]);
	if($Stringlength > 60){
		$LineTracker++;
		$printout .= "\n\t";
		$Stringlength = 1 + length($NexusTaxonSet[$khg - 1]);
		}
	$printout .=  " $NexusTaxonSet[$khg - 1]";
	}
print NEWNEXUS "$printout;\nend;";


print NEWNEXUS "\n\nBegin Characters\;\n";
print NEWNEXUS "\tDimensions nchar=$NCHAR\;";

#find which charset to use
@NexusCharSet = ();

#if transposed matrix, charnames from matrix have priority...
if($Inverted == 1){
	@NexusCharSet = @FromMatrixChars;
	}
#in non-transposed matrices:
#if there are no non-defaults, use defaults from @CharLabels...
elsif(($NonDefaultChars == 0) and ($FromMatrixCharTest == 0)){
	@NexusCharSet = @CharLabels;
	}
#if there are non-default charlabels, use them...
elsif($NonDefaultChars == 1){
	@NexusCharSet = @CharLabels;
	}
#otherwise, if there are non-default matrix labels, use them...
else{
	@NexusCharSet = @FromMatrixChars;
	}

#print charlables
unless($Inverted == 1){
	print NEWNEXUS "\n\tCharlabels\n\t";
	$LineTracker = 1;
	$printout = "";
	$Stringlength = 0;
	for($khg = 1; ($khg - 1) <= $#NexusCharSet; $khg++){
		$Stringlength = $Stringlength + 1 + length($NexusCharSet[$khg - 1]);

		if($Stringlength > 60){
			$LineTracker++;
			$printout .=  "\n\t";
			$Stringlength = 1 + length($NexusCharSet[$khg - 1]);
			}
		$printout .= " $NexusCharSet[$khg - 1]";
		}
	print NEWNEXUS "$printout;\n";
	}

#print rest of characters block
print NEWNEXUS "\n\tFormat\n\t Symbols=\" @SymbolsII \"\n\t Missing=$Missing";
if($FoundGap == 1){
	print NEWNEXUS "\n\t Gap=$Gap \[treated as $GapMode\]";
	}
if($DataType ne "standard"){
	print NEWNEXUS "\n\t DataType=$DataType ";
	}

if($Inverted == 1){
	print NEWNEXUS "\n\t Transpose";
	for($zi =1; $zi <= $NCHAR; $zi++){
		for($wq =0; $wq <= $#UninfNexusArray; $wq++){
			$SmallTemp = $UninfNexusArray[$wq][$zi];
			push(@LargeTemp, $SmallTemp);
			}
		unshift(@LargeTemp, $NexusCharSet[$zi - 1]);

		push(@InvertArray, [@LargeTemp]);
		@LargeTemp = ();
		}
	@UninfNexusArray = @InvertArray;

	$Temp = $NTAX;
	$NTAX = $NCHAR;
	$NCHAR = $Temp;
	$CharTaxCol = "Chars";
	$CharTaxRow = "Taxa ";
	}


#get max char lengths for printing
@LengthNexusArray = @UninfNexusArray;
&MaxCharLengths;

if($InterNewNexus == 1){
	print NEWNEXUS "\n\t Interleave"
	}
print NEWNEXUS "\;\n";
print NEWNEXUS "\nMatrix\n";
$BrakishO = "\[";
$BrakishC = "\]";

#Print Nexus itself:
$PrintTo = "NEWNEXUS";
$InterleafPrintLengthCheck = 0;
if($InterNewNexus == 1){
	$Explain = "\[The matrix is large; for clarity it is printed in blocks of 50 Characters\]\n";
	@NexusArrayINTER = @UninfNexusArray;
	&INTERLEAFPRINT;
	}
else{
	@NexusArrayTemp = @UninfNexusArray;
	$StartChar = 1;
	$PrintNCHAR = $NCHAR;
	&PrintNexus
	}

print NEWNEXUS "\;\nend;\n\n";

#return inverted matrix to rights
if($Inverted == 1){
	@InvertArray = ();
	@LargeTemp = ();
	for($zi =1; $zi <= $NCHAR; $zi++){
		for($wq =0; $wq <= $#UninfNexusArray; $wq++){
			$SmallTemp = $UninfNexusArray[$wq][$zi];
			push(@LargeTemp, $SmallTemp);
			}
		push(@InvertArray, [@LargeTemp]);
		@LargeTemp = ();
		}
	@UninfNexusArray = @InvertArray;
	$Temp = $NTAX;
	$NTAX = $NCHAR;
	$NCHAR = $Temp;
	}


$User_Name_Check = 1;
$User_Name_Equiv = 1;
$PerlEQ_Set = 1;

if(($AreThereTaxSets == 1) or ($AreThereCharSets == 1)){
	print NEWNEXUS "begin sets;";
	#print user taxsets
	for($rtg = 0; $rtg <= $#TaxSetFinalList; $rtg++){
		@TempTex = @{$TaxSetFinalList[$rtg]};
		$LineTracker = 1;
		$printout = "";
		$Stringlength = 0;

		if($TempTex[0] =~ /^user_deleted/){$User_Name_Check++}
		if($TempTex[0] =~ /^equivalents/){$User_Name_Equiv++}

		print NEWNEXUS "\n\ttaxset $TempTex[0] =";
		for($khg = 1; $khg <= $#TempTex; $khg++){
			if($TempTex[$khg] == 0){next}
			$Stringlength = $Stringlength + 1 + length($FinalNameList[$khg - 1]);
			if($Stringlength > 60){
				$LineTracker++;
				$printout .= "\n\t";
				$Stringlength = 1 + length($FinalNameList[$khg - 1]);
				}

			$printout .= " $FinalNameList[$khg - 1]";
			}

		if($LineTracker > 1){print NEWNEXUS "\n\t"}
		print NEWNEXUS "$printout;";
		}

	if($SomeDeleted == 1){
		$LineTracker = 1;
		$printout = "";
		$Stringlength = 0;
		print NEWNEXUS "\n\ttaxset user_deleted";
		if($User_Name_Check > 1){
			print NEWNEXUS "_$User_Name_Check";
			$Latest_Name = "user_deleted" . "_$User_Name_Check";
			}
		else{$Latest_Name = "user_deleted"}
		print NEWNEXUS " =";

		for($rtg = 0; $rtg <= $#DeletionFinalList; $rtg++){
			if($DeletionFinalList[$rtg] == 0){
				$Stringlength = $Stringlength + 1 + length($FinalNameList[$rtg]);
				if($Stringlength > 60){
					$LineTracker++;
					$printout .= "\n\t";
					$Stringlength = 1 + length($FinalNameList[$rtg]);
					}

				$printout .=  " $FinalNameList[$rtg]";
				}
			}
		if($LineTracker > 1){print NEWNEXUS "\n\t";}
		print NEWNEXUS "$printout;";
		}

	if($Is_there_anything > 0){
		$LineTracker = 1;
		$printout = "";
		$Stringlength = 0;

		print NEWNEXUS "\n\ttaxset equivalents";
		if($User_Name_Equiv > 1){
			print NEWNEXUS "_$User_Name_Equiv";
			$Latest_Equiv_Name = "equivalents" . "_$User_Name_Equiv";
			}
		else{$Latest_Equiv_Name = "equivalents"}
		print NEWNEXUS " =";

		for($RTV = 0; $RTV <= $#Tit; $RTV++){
			$NaNa = $SuggestDel{"$Tit[$RTV]"};
			if($NaNa == 0){
				$TransformedName = $Tit[$RTV];
				$TransformedName =~ s/^\<I\>//i;
				$TransformedName =~ s/\<\/I\>//i;

				$Stringlength = $Stringlength + 1 + length($TransformedName);
				if($Stringlength > 60){
					$Stringlength = 1 + length($TransformedName);
					$LineTracker++;
					$printout .= "\n\t";
					}

				$printout .= " $TransformedName";
				}
			}

		if($LineTracker > 1){print NEWNEXUS "\n\t";}
		print NEWNEXUS "$printout;";
		}

	#print user charsets
	$usex = 0;
	if(@CharsetOnlyList > 0){
		if($Totter >= 1){
			$usex = 1;
			$LineTracker = 1;
			$printout = "";
			$Stringlength = 0;
			print NEWNEXUS "\n\tcharset user_excluded =";
			for($khg = 0; $khg <= $#FinalExSet; $khg++){
				$This = $FinalExSet[$khg];
				$That = $FinalExSet[$khg + 1];
				$TheOther = $FinalExSet[$khg + 2];
				if($FinalExSet[$khg] == 1){
					$Stringlength = $Stringlength + 1 + length($khg);
					if($Stringlength >= 60){
						$LineTracker++;
						$printout .= "\n\t";
						$Stringlength = 1 + length($khg);
						}
					$printout .= " " . ($khg + 1);

					if(($That == 1) and ($TheOther == 1)){
						$printout .= "-";
						$Adder = 3;
						$That = $khg + ($Adder - 1);
						$TheOther = $FinalExSet[$khg + $Adder];
						until($TheOther != 1){
							$Adder++;
							$That = $khg + ($Adder - 1);
							$TheOther = $FinalExSet[$khg + $Adder];
							}
						$printout .= ($That + 1);
						$Stringlength = $Stringlength + 1 + length(($That + 1));
						$khg = $khg + ($Adder - 1);
						}
					}
				}
			if($LineTracker > 1){
				print NEWNEXUS "\n\t";
				};
			print NEWNEXUS "$printout\;";
			}


		for($rtg = 0; $rtg <= $#CharsetOnlyList; $rtg++){
			@TempTex = @{$CharsetOnlyList[$rtg]};
			if($TempTex[0] =~ /^PerlEQ_Uninformative/){$PerlEQ_Set++}
			if($TempTex[0] =~ /^user_excluded/){next}

			$LineTracker = 1;
			$printout = "";
			$Stringlength = 0;
			print NEWNEXUS "\n\tcharset $TempTex[0] =";
			for($khg = 1; $khg <= $#TempTex; $khg++){
				$This = $TempTex[$khg];
				$That = $TempTex[$khg + 1];
				$TheOther = $TempTex[$khg + 2];
				if($TempTex[$khg] == 1){
					$Stringlength = $Stringlength + 1 + length($khg);
					if($Stringlength >= 60){
						$LineTracker++;
						$printout .= "\n\t";
						$Stringlength = 1 + length($khg);
						}
					$printout .= " $khg";

					if(($That == 1) and ($TheOther == 1)){
						$printout .= "-";
						$Adder = 3;
						$That = $khg + ($Adder - 1);
						$TheOther = $TempTex[$khg + $Adder];
						until($TheOther != 1){
							$Adder++;
							$That = $khg + ($Adder - 1);
							$TheOther = $TempTex[$khg + $Adder];
							}
						$printout .= "$That";
						$Stringlength = $Stringlength + 1 + length($That);
						$khg = $khg + ($Adder - 1);
						}
					}
				}
			if($LineTracker > 1){
				print NEWNEXUS "\n\t";
				};
			print NEWNEXUS "$printout\;";
			}
		}

	$Latest_PerlEQ = "none";
	if($BadYes == 1){
		print NEWNEXUS "\n\tcharset PerlEQ_Uninformative";
		if($PerlEQ_Set > 1){
			print NEWNEXUS "_$PerlEQ_Set";
			$Latest_PerlEQ = "PerlEQ_Uninformative" . "_$PerlEQ_Set";
			}
		else{$Latest_PerlEQ = "PerlEQ_Uninformative"}
		print NEWNEXUS " =";
		$LineTracker = 0;
		$printout = "";
		$Stringlength = 0;
		$printadd = "";

		#prints out a compressed list, with 60chr line-length
		for($TT = 0; $TT <= $#Badlist; $TT++){
			$This = $Badlist[$TT];
			$That = $Badlist[$TT + 2];
			$Adder = 2;

			if($That == ($This + $Adder)){
				$printadd =  " $This-";
				$Adder++;
				$TheOther = $Badlist[$TT + $Adder];
				until($TheOther !=  ($This + $Adder)){
					$That = $TheOther;
					$Adder++;
					$TheOther = $Badlist[$TT + $Adder];
					}
				$printadd .= "$That";
				$TT = $TT + ($Adder - 1);
				}
			else{
				$printadd = " $This";
				}
			$Stringlength = $Stringlength + length($printadd);

			if($Stringlength >= 60){
				$printout .= "\n\t ";
				$printadd =~ s/^ //;
				$Stringlength = length($printadd);
				$LineTracker++;
				}
			$printout .= $printadd;
			}

		if($LineTracker > 1){
			print NEWNEXUS "\n\t";
			};
		print NEWNEXUS "$printout\;\n";
		}
	print NEWNEXUS "\nend;\n\n";
	}

if($AreThereCodons == 1){
	print NEWNEXUS "begin Codons;";
	@eachtype = ('n', '1', '2', '3', '?');
	$runaround = 0;
	$printout = "";
	foreach $GTO (@eachtype){
		@tempType = ();
		for($BBCITV = 0; $BBCITV <= $#UsedCodons; $BBCITV++){
			if($UsedCodons[$BBCITV] eq $GTO){
				push(@tempType, ($BBCITV + 1));
				}
			}
		unless(@tempType == 0){
			$runaround++;
			if($runaround > 1){
				$printout .= ", \n\t";
				}
			$printout .= " $GTO" . ":";
			if($#tempType == $#UsedCodons){
				$printout .=  " all";
				last;
				}

			for($TVWest = 0; $TVWest <= $#tempType; $TVWest++){
				$This = $tempType[$TVWest];
				$That = $tempType[$TVWest + 2];
				$Again  = $tempType[$TVWest + 1];

				$printout .=  " $This";
				$Adder = 2;
				$Triplet = 6;
				if($That == ($This + $Adder)){
					$printout .=  "-";
					$Adder++;

					$TheOther = $tempType[$TVWest + $Adder];
					until($TheOther !=  ($This + $Adder)){
						$That = $TheOther;
						$Adder++;
						$TheOther = $tempType[$TVWest + $Adder];
						}
					$printout .= "$That";
					$TVWest = $TVWest + ($Adder - 1);
					}

				elsif(($That == ($This + $Triplet))and ($Again == ($This + 3))){
					$printout .=  "-";
					$Adder++;
					$Triplet += 3;
					$TheOther = $tempType[$TVWest + $Adder];
					until($TheOther !=  ($This + $Triplet)){
						$That = $TheOther;
						$Adder++;
						$Triplet += 3;
						$TheOther = $tempType[$TVWest + $Adder];
						}
					$printout .= "$That\\3";
					$TVWest = $TVWest + ($Adder - 1);
					}
				}
			}
		}
	if($runaround > 1){print NEWNEXUS "\n\t";}
	print NEWNEXUS "$printout;";
	print NEWNEXUS "\nend;\n\n";
	}

print NEWNEXUS "begin PAUP;";
print NEWNEXUS "\n\texclude constant";
if($usex == 1){
	print NEWNEXUS " user_excluded";
	}
unless($Latest_PerlEQ eq "none"){
	print NEWNEXUS " $Latest_PerlEQ";
	}
print NEWNEXUS "\/only;";
print NEWNEXUS "\n\tCType";
@eachtype = ('ord', 'unord', 'dollo.up', 'dollo.dn', 'irrev.up', 'irrev.dn');
$runaround = 0;
$printout = "";
foreach $GTO (@eachtype){

	@tempType = ();
	for($BBCITV = 0; $BBCITV <= $#CharTypes2; $BBCITV++){
		if($CharTypes2[$BBCITV] eq $GTO){
			push(@tempType, ($BBCITV + 1));
			}
		}
	unless(@tempType == 0){
		$runaround++;
		if($runaround > 1){
			$printout .= ", \n\t";
			}
		$printout .= " $GTO" . ":";
		if($#tempType == $#CharTypes2){
			$printout .=  " all";
			last;
			}

		for($TVWest = 0; $TVWest <= $#tempType; $TVWest++){
			$This = $tempType[$TVWest];
			$That = $tempType[$TVWest + 2];

			$printout .=  " $This";
			$Adder = 2;
			if($That == ($This + $Adder)){
				$printout .=  "-";
				$Adder++;

				$TheOther = $tempType[$TVWest + $Adder];
				until($TheOther !=  ($This + $Adder)){
					$That = $TheOther;
					$Adder++;
					$TheOther = $tempType[$TVWest + $Adder];
					}
				$printout .= "$That";
				$TVWest = $TVWest + ($Adder - 1);
				}
			}

		}
	}
if($runaround > 1){print NEWNEXUS "\n\t";}
print NEWNEXUS "$printout;";


if(($SomeDeleted == 1) or ($Is_there_anything > 0) or ($AreTherePolym == 1)){
	if(($SomeDeleted == 1) or ($Is_there_anything > 0)){
		print NEWNEXUS "\n\tdelete";
		if($SomeDeleted == 1){
			print NEWNEXUS " $Latest_Name";
			}
		if($Is_there_anything > 0){
			print NEWNEXUS " $Latest_Equiv_Name";
			}
		print NEWNEXUS "\/only;";
		}

	if($AreTherePolym == 1){
		print NEWNEXUS "\n\tPSet MSTaxa=$MStaxa;";
		}
	}
print NEWNEXUS "\nend;\n\n";


print NEWNEXUS "begin assumptions\;\n";
if($FoundGap == 1){
	print NEWNEXUS "\tOptions GapMode = $GapMode;\n";
	}

print NEWNEXUS "\ttypeset ";

if($UsedTypeset !~ m/^\*/){print NEWNEXUS "\*"}
print NEWNEXUS "$UsedTypeset =";
$runaround = 0;
$printout = "";
foreach $GTO (@eachtype){

	@tempType = ();
	for($BBCITV = 0; $BBCITV <= $#CharTypes2; $BBCITV++){
		if($CharTypes2[$BBCITV] eq $GTO){
			push(@tempType, ($BBCITV + 1));
			}
		}
	unless(@tempType == 0){
		$runaround++;
		if($runaround > 1){
			$printout .= ", \n\t";
			}
		$printout .= " $GTO" . ":";
		if($#tempType == $#CharTypes2){
			$printout .=  " all";
			last;
			}

		for($TVWest = 0; $TVWest <= $#tempType; $TVWest++){
			$This = $tempType[$TVWest];
			$That = $tempType[$TVWest + 2];

			$printout .=  " $This";
			$Adder = 2;
			if($That == ($This + $Adder)){
				$printout .=  "-";
				$Adder++;

				$TheOther = $tempType[$TVWest + $Adder];
				until($TheOther !=  ($This + $Adder)){
					$That = $TheOther;
					$Adder++;
					$TheOther = $tempType[$TVWest + $Adder];
					}
				$printout .= "$That";
				$TVWest = $TVWest + ($Adder - 1);
				}
			}

		}

	}
if($runaround > 1){print NEWNEXUS "\n\t";}

print NEWNEXUS "$printout;\nend;";

#Program terminates here
#below is the run-time calculation
@EndTime = localtime(time);

if($EndTime[0] -  $StartTime[0] >= 0){
	$RunSecs = $EndTime[0] -  $StartTime[0]
	}
else{
	$RunSecs = ($EndTime[0] -  $StartTime[0]) + 60;
	$EndTime[1]--;
	}

if($EndTime[1] -  $StartTime[1] >= 0){
	$RunMins = $EndTime[1] -  $StartTime[1]
	}
else{
	$RunMins = ($EndTime[1] -  $StartTime[1]) + 60;
	$EndTime[2]--;
	}

if($EndTime[2] -  $StartTime[2] >= 0){
	$RunHours = $EndTime[2] -  $StartTime[2]
	}
else{
	$RunHours = ($EndTime[2] -  $StartTime[2]) + 24;
	$EndTime[7]--;
	}

if($EndTime[7] -  $StartTime[7] >= 0){
	$RunDays = $EndTime[7] -  $StartTime[7]
	}
else{
	$RunDays = ($EndTime[7] -  $StartTime[7]) + 365
	}

print "\n\nRun took $RunDays" . "d, $RunHours" . "h, $RunMins" . "m, $RunSecs" . "s\n";

print OUTPUT "<br><br><hr>\n";
print OUTPUT "\t\t<b>Run took $RunDays" . "d, $RunHours" . "h, $RunMins" . "m, $RunSecs" . "s";
print OUTPUT "\n</b><br><hr></p>\n";
print OUTPUT "\n<\/body>\n<\/html>";

print "\nHope it worked!\n\n";





#############################################
#############################################
###					  ###
###	 SSSS	UU  UU	BBBBB	 SSSS	  ###
###	SSSSSS	UU  UU	BBBBBB	SSSSSS	  ###
###	SS  SS	UU  UU	BB  BB	SS  SS	  ###
###	 SS	UU  UU	BBBBB	 SS	  ###
###	  SS	UU  UU	BBBBBB	  SS	  ###
###	   SS	UU  UU	BB   BB	   SS	  ###
###	SS  SS	UU  UU	BB   BB	SS  SS	  ###
###	SSSSSS	UUUUUU	BBBBBBB	SSSSSS	  ###
###	 SSSS	 UUUU	BBBBBB	 SSSS	  ###
###					  ###
#############################################
#############################################


#########################################################
#Sub to print time
sub Timeline {
	my @T2 = localtime(time);

	if($T2[0] -  $T1[0] >= 0){
		$TSecs = $T2[0] -  $T1[0]
		}
	else{
		$TSecs = ($T2[0] -  $T1[0]) + 60;
		$T2[1]--;
		}

	if($T2[1] -  $T1[1] >= 0){
		$TMins = $T2[1] -  $T1[1]
		}
	else{
		$TMins = ($T2[1] -  $T1[1]) + 60;
		$T2[2]--;
		}

	if($T2[2] -  $T1[2] >= 0){
		$THours = $T2[2] -  $T1[2];
		}
	else{
		$THours = ($T2[2] -  $T1[2]) + 24;
		}
	my @Output = ($TSecs, $TMins, $THours);
	foreach my $x (@Output){
		if(length($x) < 2){
			until(length($x) == 2){
				$x = '0' . $x;
				}
			}
		}

	my $Percentum = shift;
	if(length($Percentum) < 3){
		until(length($Percentum) == 3){
			$Percentum = " " . $Percentum;
			}
		}
	$Percentum = $Percentum  . "\%";
	push(@Output, $Percentum);
	@Output;
	}
#end of sub to print time
#########################################################


#########################################################
#Sub to recode ordered multistate chars so that they have
#all intermediate states...
sub OrderRecode {

	for($XX = 1; $XX <= $NCHAR; $XX++){
		$TypeofChar = $CharTypes2[$XX - 1];
		for($YY = 0; $YY < $NTAX; $YY++){
			$Coding = $RecodeArray[$YY][$XX];
			if(($TypeofChar ne "unord") and ($Coding =~ /\)$|\}$/)){
				@Temp = split(//, $Coding);
				$First = $Temp[1];
				$Last = $Temp[$#Temp - 1];
				$Middle = "";
				$Startup = 0;
				foreach $ZZ (@RecodeSymbols){
					if($ZZ eq $First){
						$Startup = 1}
					if($Startup == 1){
						$Middle .= $ZZ;
						}
					if($ZZ eq $Last){
						$Startup = 0}
					}
				if($Coding =~ /\}$/){
					$Middle = "{" . $Middle . "}";
					}
				else{
					$Middle = "(" . $Middle . ")";
					}
				$RecodeArray[$YY][$XX] = $Middle;
				}
			}
		}
	}

#end of sub to recode ordered multistate chars
#########################################################


#########################################################
#Sub to assess the overall informativeness of a character
sub OverallInform {
	for($f = 0; $f < $NCHAR; $f++){
		$TypeofChar = $CharTypes2[$f];
		@TempStateData = @{$CharStatesInf[$f]};

		#get rid of 'missing' and gap (if gapmode is 'missing')
		pop(@TempStateData);
		if(($FoundGap == 1) and ($GapMode eq "missing")){
			pop(@TempStateData);
			}

		$BadCheck = 0;
		@KillerList = ();
		&Redistrib;

		if(@KillerList >= 2){
			$BadCheck = 2;
			}
		else{
			unless($FinalExSet[$f] == 1){&BadChars;}
			}
		}
	}
#end of sub to assess the overall informativeness
#########################################################


#########################################################
#little sub judge if overall score data is informative
sub Redistrib {
	#find scores with 2 or more codings
	@KillerList = ();

	if($TypeofChar =~ /unord/){
		for($UY = 0; $UY <= $#TempStateData; $UY++){
			if($TempStateData[$UY] >= 2){
				push(@KillerList, $InfSymbols2[$UY]);
				}
			}
		}

	if($TypeofChar =~ /^ord$|dollo|irrev\.up/){
		for($UY = $#TempStateData; $UY >= 0; $UY--){
			if($TempStateData[$UY] >= 2){
				push(@KillerList, $InfSymbols2[$UY]);
				last;
				}
			}
		}

	if($TypeofChar =~ /^ord$|dollo|irrev\.dn/){
		for($UY = 0; $UY <= $#TempStateData; $UY++){
			if($TempStateData[$UY] >= 2){
				push(@KillerList, $InfSymbols2[$UY]);
				last;
				}
			}
		}

	#for irrevesable chars add base-state
	if($TypeofChar =~ /irrev/){
		if($TypeofChar =~ /\.up/){
			for($UY = 0; $UY <= $#TempStateData; $UY++){
				if($TempStateData[$UY] > 0){
					push(@KillerList, $InfSymbols2[$UY]);
					last;
					}
				}
			}
		if($TypeofChar =~ /\.dn/){
			for($UY = $#TempStateData; $UY >= 0; $UY--){
				if($TempStateData[$UY] > 0){
					push(@KillerList, $InfSymbols2[$UY]);
					last;
					}
				}
			}
		}

	#make sure the same char isn't supporting the upper and lower bounds
	@New = ();
	$Tempestu = $f + 1;
	foreach $i (@KillerList){

		unless(grep(/$i/, @New) == 1){
			push(@New, $i);
			}
		}
	@KillerList = @New;
	}

#end of little sub to redistriute poly scores
#########################################################


#########################################################
#Sub to replace uninformative unordered data
sub ReplaceUnord {
	@Replacelist = ();
	@ReplacelistTo = ();
	$GoReplace = 0;

	#get state data for the char in question
	@TempSymbols2 = @SymbolsII;
	if(($FoundGap == 1) and ($GapMode ne "missing")){
		push(@TempSymbols2, $Gap);
		}

	#first bring down any polymorph to a single code if poss.

	for($z = 0; $z <= $#TempSymbols2; $z++){
		$pushit = 0;
		if($TempStateData[$z] == 0){
			foreach $parrot (@TempPolyData){
				if($parrot =~ /$TempSymbols2[$z]/){
					$tempin = $TempSymbols2[$z];
					$pushit = 1;
					$DoAgain = 1;
					}
				}
			}
		if($pushit == 1){
			push(@Replacelist, $tempin);
			push(@ReplacelistTo, $Missing);
			}
		}

	unless(@Replacelist > 0){
		#for each state, check it has 2 or more codings
		for($z = 0; $z <= $#TempSymbols2; $z++){
			$pushit = 0;

			if($TempStateData[$z] == 1){
				$tempin = $TempSymbols2[$z];
				$pushit = 1;
				$DoAgain = 1;
				}

			unless($pushit == 0){
				push(@Replacelist, $tempin);
				push(@ReplacelistTo, $Missing);
				$DoAgain = 1;
				}
			}
		}

	if($DoAgain == 1){
		$CurrentChar = $i + 1;
		&Replacing;
		}
	}
#end of sub to replace uninformative unordered data
#########################################################


#########################################################
#Sub Lower bounds
sub ReplaceLower {
	@Replacelist = ();
	@ReplacelistTo = ();
	$GoReplace = 0;

	#get state data for the char in question
	@TempSymbols2 = @SymbolsII;
	if(($FoundGap == 1) and ($GapMode ne "missing")){
		push(@TempSymbols2, $Gap);
		}

	#check lower bounds
	$LoopStop = 0;

	#find if lowest state 2 or more codings
	for($z = 0; $z <= $#TempSymbols2; $z++){
		$pushit = 0;

		if($LoopStop == 1){last}

		#if no regular coding, delete from polym.
		elsif($TempStateData[$z] == 0){
			foreach $parrot (@TempPolyData){
				if($parrot =~ /$TempSymbols2[$z]/){
					$tempin = $TempSymbols2[$z];
					if($z == $#TempSymbols2){$Nexttempin = $Missing}
					else{$Nexttempin = $TempSymbols2[$z + 1];}
					$pushit = 1;
					$DoAgain = 1;
					}
				}
			}

		#if one regular coding, delete.
		elsif($TempStateData[$z] == 1){
			$pushit = 1;
			$tempin = $TempSymbols2[$z];
			if($z == $#TempSymbols2){$Nexttempin = $Missing}
			else{$Nexttempin = $TempSymbols2[$z + 1];}
			}

		#state has informative no of codings
		else{
			$LoopStop = 1;
			$pushit = 0;
			}

		unless($pushit == 0){
			push(@Replacelist, $tempin);
			push(@ReplacelistTo, $Nexttempin);
			$LoopStop = 1;
			$DoAgain = 1;
			}
		}
	if($DoAgain == 1){
		$CurrentChar = $i + 1;
		&Replacing
		}
	}
#end of Sub Lower bounds
#########################################################


#########################################################
#Sub Upper bounds
sub ReplaceUpper {
	@Replacelist = ();
	@ReplacelistTo = ();
	$GoReplace = 0;

	#get state data for the char in question
	@TempSymbols2 = @SymbolsII;
	if(($FoundGap == 1) and ($GapMode ne "missing")){
		push(@TempSymbols2, $Gap);
		}

	#check upper bounds
	$LoopStop = 0;

	#find if highest state 2 or more codings
	for($z = $#TempSymbols2; $z >= 0; $z--){
		$pushit = 0;

		if($LoopStop == 1){
			last;
			}

		#if no regular coding, delete from polym.
		elsif($TempStateData[$z] == 0){
			foreach $parrot (@TempPolyData){
				if($parrot =~ /$TempSymbols2[$z]/){
					$tempin = $TempSymbols2[$z];
					if($z == 0){$Nexttempin = $Missing}
					else{$Nexttempin = $TempSymbols2[$z - 1];}
					$pushit = 1;
					$DoAgain = 1
					}
				}
			}

		#if one regular coding, delete.
		elsif($TempStateData[$z] == 1){
			$pushit = 1;
			$tempin = $TempSymbols2[$z];
			if($z == 0){$Nexttempin = $Missing}
			else{$Nexttempin = $TempSymbols2[$z - 1];}
			}
		else{
			$LoopStop = 1;
			$pushit = 0;
			}

		unless($pushit == 0){

			push(@Replacelist, $tempin);
			push(@ReplacelistTo, $Nexttempin);
			$LoopStop = 1;
			$DoAgain = 1
			}
		}

	if($DoAgain == 1){
		$CurrentChar = $i + 1;
		&Replacing
		}
	}
#end of Sub Upper bounds
#########################################################


#########################################################
#Sub to replace from the Replacelist
sub Replacing {
	$HasAnyThingChanged = 1;
	for($p = 0; $p <= $#Replacelist; $p++){
		$CharFrom = $Replacelist[$p];
		$CharTo = $ReplacelistTo[$p];

		for($b = 0; $b < $NTAX; $b++){
			$TestState = $UninfNexusArray[$b][$CurrentChar];
			if(($TestState eq $CharFrom) or ($TestState =~ /\)$|\}$/)){

				#if a regular char
				if($TestState eq $CharFrom){
					$UninfNexusArray[$b][$CurrentChar] = $CharTo;
					}

				#if polymorph
				elsif($TestState =~ /\)$|\}$/){

					#if replacement is "?" delete from brackets
					if($CharTo eq $Missing){
						$UninfNexusArray[$b][$CurrentChar] =~ s/$CharFrom//;
						}
					#otherwise replace
					else{
						$UninfNexusArray[$b][$CurrentChar] =~ s/$CharFrom/$CharTo/;

						#make sure that the polym doesn't have the same state twice
						@SplitList = split(//, $UninfNexusArray[$b][$CurrentChar]);
						$FinalString = "";
						foreach $Additional (@SplitList){
							#rebuild polym char without doubles
							unless($FinalString =~ quotemeta($Additional)){
								$FinalString .= $Additional;
								}
							}
						$UninfNexusArray[$b][$CurrentChar] = $FinalString;
						}

					#if only brackets left, replace with "?"
					if(length($UninfNexusArray[$b][$CurrentChar]) == 2){
						$UninfNexusArray[$b][$CurrentChar] = $Missing;
						}
					#if only one char replace with the single char
					if(length($UninfNexusArray[$b][$CurrentChar]) == 3){
						$UninfNexusArray[$b][$CurrentChar] =~ s/\{|\(//;
						$UninfNexusArray[$b][$CurrentChar] =~ s/\}|\)//;
						}
					}

				#check to see if that char has already been changed in a previous round
				for($rD = 0; $rD <= $#ChangedData; $rD++){
					@tvu = @{$ChangedData[$rD]};
					$TipTest = 0;
					if(
						($tvu[0] == $b) and
						($tvu[1] == $CurrentChar) and
						($tvu[3] == $TestState)
						){
						$TipTest = 1;
						$ChangedData[$rD] = [($b, $CurrentChar, $tvu[2], $CharTo)];
						last;
						}
					}
				unless($TipTest == 1){
					@TempList = ($b, $CurrentChar, $TestState, $CharTo);
					push(@ChangedData, [@TempList]);
					}

				$TestChangedData = 1;
				}
			}
		}
	}
#end of sub to replace from the Replacelist
#########################################################


#########################################################
#little sub to sort exclusion list the order in the original NEXUS
sub TitSorter {
	for($FBS = 0; $FBS <= $#NexusArray; $FBS++){
		if($NexusArray[$FBS][0] eq $a){$This_n = $FBS}
		if($NexusArray[$FBS][0] eq $b){$That_n = $FBS}
		}
	return($This_n <=> $That_n);
	}
#end of sub
#########################################################


#########################################################
#Sub to get print changes
sub PrintChange {
	@TempChanged = @{$ChangedData[$g]};
	$SpeciesTemp =  $UninfNexusArray[$TempChanged[0]][0];
	$SpeciesTemp =~ s/^\'//;
	$SpeciesTemp =~ s/\'$//;
	$ChTemp =  $TempChanged[1];

	if($TempChanged[2] eq $Gap){
		$StateTemp = "Gap symbol \'$Gap\'";
		}
	else{$StateTemp = "State \'$TempChanged[2]\'";}

	print OUTPUT "\n<tr>\n\t<td class=\"nexusSp\"><i>$SpeciesTemp</i></td>";
	print OUTPUT "\n\t<td class=\"nexusSp\"><div align=\"center\">$ChTemp</div></td>";
	print OUTPUT "\n\t<td class=\"nexusSp\">&nbsp;$StateTemp uninformative, changed to \'$TempChanged[3]\'</td>\n</tr>";
	}
#end of sub to get print changes
#########################################################


#########################################################
#sub to replace character states of uninf chars
sub BadChars {
	for($i = 0; $i < $NTAX; $i++){
		$UninfNexusArray[$i][$f + 1] = $Missing;
		}
	push(@Badlist, ($f + 1)) ;
	$BadYes = 1;
	$TestChangedData++;
	}
#end of sub to replace character states of uninf chars
#########################################################


#########################################################
#Sub to print EQ data
sub PrintEQ {
	$EQLength = $NTAX - 1;
	$CommentOrNot = 0;
	%SuggestDel = ();

	# @FinalEQdata [(Indextaxon, NoSpaces, Missing, Eq1, Eq2, Eq3...),...]
	# @DeletionFinalList (1/0)
	for($DD = 0; $DD <= $EQLength; $DD++){
		$NameTempBit = $FinalEQdata[$DD][0];
		$SuggestDel{"$NameTempBit"} = 1;
		@Tit = keys(%SuggestDel);
		}

	for($DD = 0; $DD <= $EQLength; $DD++){
		$BGColor = "";

		if($DD%2 == 0){
			$BGColor = "class=\"clear\"";
			}
		@TempEQ = @{$FinalEQdata[$DD]};
		$endlist = $#TempEQ;
		$TestName = $TempEQ[0];
		$Additional = $TempEQ[1];
		$ExtendPrint = 0;
		$Dull1 = "";
		$Dulloff = "";
		if($DeletionFinalList[$DD] == 0){
			$Dull1 = " <p class=\"Grey\">";
			$Dulloff = " </p>";
			$HowManyMissing = '-'
			}
		elsif($TempEQ[2] == 0){$HowManyMissing = "none"}
		else{
			$HowManyMissing = $TempEQ[2];
			$ExtendPrint = 1;
			}
		$quoteless = $TestName;
		$quoteless =~ s/^\'//;
		$quoteless =~ s/\'$//;

		print OUTPUT "\n<tr $BGColor>\n\t<td class=\"nexusEq\">$Dull1<i>$quoteless</i>$Dulloff\</td>";
		print OUTPUT "\n\t<td class=\"nexusEq\"><div align=\"center\">$Dull1\&nbsp;\[$HowManyMissing";
		if($ExtendPrint == 1){
			$TopHat = "";
			$Spresso = "s";
			if($HowManyMissing == ($NCHAR - $EXCLUDEDcountTotal)){
				$TopHat = "<b>\^</b>";
				$CommentOrNot = 1;
				}
			if($HowManyMissing == 1){
				$Spresso = ""
				}
			print OUTPUT " char$Spresso$TopHat, ";
			$PercentMissing = ($HowManyMissing/($NCHAR - $EXCLUDEDcountTotal)) * 100;
			printf OUTPUT ("%.1f\%" , $PercentMissing);
			#if($HowManyMissing == 1){
			#	print OUTPUT "<what's this for?>";
			#	}
			}
		print OUTPUT "\]&nbsp;$Dulloff</div></td>";
		if($DeletionFinalList[$DD] == 0){
			$printout = "\t<td class=\"nexusEq\"><p class=\"Grey\">Taxon deleted by user</p>";
			$SuggestDel{"$TestName"} = 2;
			}
		else{
			$LineTracker = 1;
			$printout = "\t<td class=\"nexusEq\">";
			$Stringlength = 0;

			for($AA = 3; $AA <= $endlist; $AA++){
				$EQstuff = $TempEQ[$AA];
				$quoteless = $TempEQ[$AA];
				$quoteless =~ s/^\<i\>\'/\<i\>/i;
				$quoteless =~ s/\'\<\/i\>/\<\/i\>/i;



				if($AA < $endlist){
					$quoteless .= ", ";
					}

				if($Stringlength + (length($quoteless)) >= 70){
					$spacer = "</td>\n</tr>\n<tr $BGColor>\n\t<td class=\"nexusEq\">&nbsp;</td>\n\t<td class=\"nexusEq\">&nbsp;</td>\n\t<td class=\"nexusEq\">";
					$printout .= $spacer;
					$Stringlength = length($quoteless);
					}
				else{$Stringlength = $Stringlength + length($quoteless);}
				$printout .= $quoteless;

				if($SuggestDel{"$TestName"} == 1){
					&Delete_List_Sorter;
					}
				}
			}
		print OUTPUT "\n$printout</td>\n</tr>";
		}
	print OUTPUT "\n</table><p><br>";
	if($CommentOrNot == 1){
		print OUTPUT "\n<b>\^</b>=Taxon has missing entries only<br>"
		}
	}
#End of sub to print EQ data
#########################################################


#########################################################
#Sub to create a suggested deletion list
sub Delete_List_Sorter {
	$deleteName = $EQstuff;
	$deleteName =~ s/\([A-E]\S?\)(\,? ?)$//;
	$deleteName =~ s/\<I\>//i;
	$deleteName =~ s/\<\/I\> //i;

	@Test_for_delete = split(/\(/, $EQstuff);

	if($Test_for_delete[$#Test_for_delete] =~ /A/){
		$SuggestDel{"$deleteName"} = 0;
		}
	if($Test_for_delete[$#Test_for_delete] =~ /B/){
		$SuggestDel{"$deleteName"} = 0;
		}
	if($Test_for_delete[$#Test_for_delete] =~ /C/){
		$SuggestDel{"$deleteName"} = 0;
		}
	}
#end of Sub to create a suggested deletion list
#########################################################


#########################################################
#subroutine to Compare missing data
sub Comparisons {
	$TestTotalComp = 0;
	@FinalEQdata = ();

	#1. check to see if there are any codings on that species
	for($i = 0; $i < $NTAX; $i++){
		@TestSpecies = @{$MissingNexusArray[$i]};
		$InformList[$i] = "!";
		for($g = 1; $g <= $NCHAR; $g++){
			if(($FoundGap == 1) and ($GapMode eq "missing")){
				if(ord($TestSpecies[$g] ne $Missing) and ($TestSpecies[$g] ne $Gap)){
					$InformList[$i] = "*"
					}
				}

			#if gap is a newstate
			else{
				if($TestSpecies[$g] ne $Missing){
					$InformList[$i] = "*"
					}
				}
			}
		}

	#2. Big equivalence loop
	for($i = 0; $i < $NTAX; $i++){
		if($DeletionFinalList[$i] == 0){
			@TempEQ = ();
			$species = $MissingNexusArray[$i][0];
			push(@TempEQ, $species);
			$Additional = $Namelength - length($species);
			push(@TempEQ, $Additional);
			push(@TempEQ, $NCHAR);
			push(@TempEQ, "deleted");
			push(@FinalEQdata, [@TempEQ]);
			@TempEQ = ();
			next;
			}

		#get test species data
		@TestSpecies = @{$MissingNexusArray[$i]};
		$TestName = $MissingNexusArray[$i][0];
		@TestOutput = ();
		for($x = 0; $x < $NTAX; $x++){
			if($x == $i){next}
			elsif($DeletionFinalList[$x] == 0){next}
			@CompareSpecies = @{$MissingNexusArray[$x]};
			$CompareName = $MissingNexusArray[$x][0];
			$NonEquiv = 0;
			$yMissing = 0;
			$xMissing = 0;
			$Symmetric = 1;
			$xMiss_yCode = 0;
			$xCode_yMiss = 0;
			$HowManyMissing = 0;
			@Compy = ();

			#get compare species data
			for($w = 1; $w <= $NCHAR; $w++){
				if($FinalExSet[$w - 1] == 1){
					next
					}

				$TestRaw = $TestSpecies[$w];
				$CompRaw = $CompareSpecies[$w];

				#adjust if gap = missing
				if(($FoundGap == 1) and ($GapMode eq "missing")){
					if($TestRaw eq $Gap){
						$TestRaw = $Missing;
						}
					if($CompRaw eq $Gap){
						$CompRaw = $Missing;
						}
					}

				#test polymorph code
				if($TestRaw =~ /\)$/){

					#comp a regular code
					if($CompRaw !~ /\}$|\)$/){
						if($TestRaw =~ quotemeta($CompRaw)){
							#real codes match
							$TestRaw = $CompRaw;
							}
						#otherwise they are different
						else{
							$TestRaw = 1;
							$CompRaw = 2;
							}
						}

					#comp an MS
					elsif($CompRaw =~ /\}$|\)$/){
						@Compy = split(//, $CompRaw);
						pop(@Compy);
						shift(@Compy);

						$CompMatchup = 0;
						foreach $Yu (@Compy){

							#if the sets of codes intersect...
							if($TestRaw =~ /$Yu/){
								#comp uncert - it's uncertain, but overlapping
								if($CompRaw =~ /\}$/){

									$CompRaw = $Missing;
									$CompMatchup = 1;
									last;
									}
								#comp a variable - real codes overlap
								else{
									$CompRaw = $Yu;
									$TestRaw = $Yu;
									$CompMatchup = 1;
									last;
									}
								}
							#otherwise the sets don't intersect, and the codes are different
							}
						if($CompMatchup == 0){
							$CompRaw = 1;
							$TestRaw = 2;
							}

						}
					}

				#test uncertain code
				elsif($TestRaw =~ /\}$/){
					# add the fact that the code is a form of 'missing'
					# - this avoids funny accounting when I recode the MS scores below
					$HowManyMissing++;

					#comp a regular code
					if($CompRaw !~ /\)$|\}$/){
						if(($TestRaw =~ quotemeta($CompRaw)) or ($CompRaw eq $Missing)){
							#uncertain, but overlaps real code or missing
							$TestRaw = $Missing;
							#but not really a full missing entry...
							$HowManyMissing--;
							}
						#otherwise they are different
						else{
							$TestRaw = 1;
							$CompRaw = 2;
							}


						}

					#comp a polymorphic
					elsif($CompRaw =~ /\)$|\}$/){

						@Compy = split(//, $CompRaw);
						pop(@Compy);
						shift(@Compy);


						$CompMatchup = 0;
						foreach $Yu (@Compy){
							#if sets intersect
							if($TestRaw =~ /$Yu/){

								#comp uncert - both uncertain, but overlapping
								if($CompRaw =~ /\}$/){
									$CompMatchup = 1;
									$CompRaw = $Missing;
									$TestRaw = $Missing;
									#but not really a full missing entry...
									$HowManyMissing--;

									last;
									}
								#comp a polmmorph - uncert overlaps real code
								else{
									$CompMatchup = 1;
									$TestRaw = $Missing;
									#but not really a full missing entry...
									$HowManyMissing--;

									last;
									}
								}
							#otherwise they are different
							}
						if($CompMatchup == 0){
							$TestRaw = 1;
							$CompRaw = 2;

								}

						}
					}

				#test regular code
				elsif($TestRaw !~ /\)$|\}$/){
					#if comp is MS
					if($CompRaw =~ /\)$|\}$/){
						#if codes intersect...
						if($CompRaw =~ quotemeta($TestRaw)){
							#if comp is polymorphic
							if($CompRaw =~ /\}$/){
								$CompRaw = $Missing;
								}
							else{
								$CompRaw = $TestRaw;
								}
							}

						elsif(($TestRaw eq $Missing) and ($CompRaw =~ /\}$/)){
							$CompRaw = $Missing;
							}

						elsif(($TestRaw eq $Missing) and ($CompRaw =~ /\)$/)){
							$CompRaw = 1;
							}

						else{
							$TestRaw = 1;
							$CompRaw = 2;
							}
						}
					}

				#get ASCII values for the scores
				$TestOrd = ord($TestRaw);
				$CompOrd = ord($CompRaw);
				$TestTotalComp++;

				#comparisons
				if(($TestOrd != ord($Missing)) and ($CompOrd != ord($Missing))){
					if($TestOrd != $CompOrd){
						$NonEquiv = 1;
						}
					}
				elsif(($TestOrd == ord($Missing)) and ($CompOrd != ord($Missing))){
					$TotalInvalid++;
					$HowManyMissing++;
					$xMissing = 1;
					$xMiss_yCode = 1;
					$Symmetric = 0;
					}
				elsif(($TestOrd != ord($Missing)) and ($CompOrd == ord($Missing))){
					$TotalInvalid++;
					$yMissing = 1;
					$xCode_yMiss = 1;
					$Symmetric = 0;

					}
				elsif(($TestOrd == ord($Missing)) and ($CompOrd == ord($Missing))){
					$HowManyMissing++;
					$TotalInvalid++;
					$yMissing = 1;
					$xMissing = 1;

					}
				}

			if($NonEquiv == 1){next}
			elsif($NonEquiv == 0){
				if($Symmetric == 1){
					if(($xMissing == 0) and ($yMissing == 0)){
						push(@TestOutput, "$CompareName(A)");
						next
						}
					else{
						push(@TestOutput, "$CompareName(B)");
						next
						}
					}
				elsif($Symmetric == 0){
					if(($xMissing == 0) and ($yMissing == 1)){
						push(@TestOutput, "$CompareName(C)");
						next
						}
					elsif(($xMissing == 1) and ($yMissing == 0)){
						push(@TestOutput, "$CompareName(E)");
						next
						}
					elsif(($xMissing == 1) and ($yMissing == 1)){
						if(($xCode_yMiss == 1) and ($xMiss_yCode == 1)){
							push(@TestOutput, "$CompareName(D)");
							next
							}
						elsif($xMiss_yCode == 0){
							push(@TestOutput, "$CompareName(C)");
							next
							}
						elsif($xCode_yMiss == 0){
							push(@TestOutput, "$CompareName(E)");
							next
							}
						}
					}
				}#end of if $NonEquiv = 0
			}#end of loop to find taxa equivalents

		@TempEQ = ();
		push(@TempEQ, $TestName);
		$Additional = $Namelength - length($TestName);
		push(@TempEQ, $Additional);
		push(@TempEQ, $HowManyMissing);
		if(@TestOutput > 0){
			push(@TempEQ, @TestOutput)
			}
		else{push(@TempEQ, "Shows no equivalence")}
		push(@FinalEQdata, [@TempEQ]);

		#bit to print-out percent progress...
		$TotUp = $Excess + int(100*(($i+1)/(2*($NTAX))));
		if($TotUp >= ($OldTotUp + 5)){
			$TimeIt+=5;
			$OldTotUp = $TotUp;
			if($OldTotUp%5 != 0){
				until($OldTotUp %5 == 0){
					$OldTotUp--;
					}
				}
			@TimePrint = Timeline($TimeIt);
			print "\n\*\* $TimePrint[3] - $TimePrint[2]h, $TimePrint[1]m, $TimePrint[0]s elapsed since last check  \*\*";
			@T1  = localtime(time);
			}
		}#end of first round of testing
	$Excess = $TotUp;
	}
#end of sub to Compare missing data
#########################################################


#########################################################
#sub to add correct suffixes to the data
sub SuffixFixer {
	# @FinalEQdata [(Indextaxon, NoSpaces, Missing, Eq1, Eq2, Eq3...),...]
	for($X = 0; $X < $NTAX; $X++){
		$NoStar = 0;
		@TestSpecies = @{$FinalEQdata[$X]};
		$TestName = $TestSpecies[0];
		@TempEquivs =();
		for($Y = 3; $Y <= $#TestSpecies; $Y++){
			push(@TempEquivs,  $TestSpecies[$Y]);
			if($TestSpecies[$Y] =~ /\(E\)$/){
				$NoStar = 1;
				}
			}

		#check equives
		for($Y = 3; $Y <= $#TestSpecies; $Y++){

			#if an A
			if($TestSpecies[$Y] =~ /\(A(\*?)\)$/){
				$FinalEQdata[$X][$Y] = "<i>" . $FinalEQdata[$X][$Y];
				$FinalEQdata[$X][$Y]=~ s/\(A\)$/<\/i> \(A\*\)/;
				}

			#if a B
			elsif($TestSpecies[$Y] =~ /\(B(\*?|\!?)\)$/){
				$Additon = "";
				#get number of equiv species
				$a = $TestSpecies[$Y];
				$a =~ s/\(B(\*?|\!?)\)$//;
				$b = "";
				&TitSorter;
				@EquivEquiv = @{$FinalEQdata[$This_n]};

				#if equiv has no data...
				if($EquivEquiv[2] == $NCHAR){
					$Additon = "\!";
					}
				#else assign star or not...
				else{
					unless($NoStar == 1){
						$Additon = "\*";
						}
					}
				$FinalEQdata[$X][$Y] = "<i>" . $FinalEQdata[$X][$Y];
				$FinalEQdata[$X][$Y] =~ s/\(B(\*?|\!?)\)$/<\/i> \(B$Additon\)/;
				}

			#if a C
			elsif($TestSpecies[$Y] =~ /\(C(\*?|\!?)\)$/){
				$Additon = "";
				#get number of equiv species
				$a = $TestSpecies[$Y];
				$a =~ s/\(C(\*?|\!?)\)$//;
				$b = "";
				&TitSorter;
				@EquivEquiv = @{$FinalEQdata[$This_n]};

				#if equiv has no data...
				if($EquivEquiv[2] == ($NCHAR - $EXCLUDEDcountTotal)){
					$Additon = "\!";
					}

				#else assign star or not...
				else{
					unless($NoStar == 1){
						$Additon = "\*";
						for($M = 3; $M <= $#EquivEquiv; $M++){
							unless($EquivEquiv[$M] =~ /\(D\S?\)$/){
								$EqivName = $EquivEquiv[$M];
								$EqivName =~ s/\([A-E]\S?\)$//;
								$Metequiv = quotemeta($EqivName);
								if((grep(/^$Metequiv\(/, @TempEquivs) == 0) and ($EqivName ne $TestName)){
									$Additon = "";
									}
								}
							}
						}
					}
				$FinalEQdata[$X][$Y] = "<i>" . $FinalEQdata[$X][$Y];
				$FinalEQdata[$X][$Y] =~ s/\(C(\*?|\!?)\)$/<\/i> \(C$Additon\)/;
				}
			else{
				unless($TestSpecies[$Y] =~ /Shows no equivalence/){
					$FinalEQdata[$X][$Y] = "<i>" . $FinalEQdata[$X][$Y];
					$FinalEQdata[$X][$Y] =~ s/\(([A-E])(\*?|\!?)\)$/<\/i> \($1\)/;
					}
				next}
			}
		}
	}
#end of sub to add correct suffixes to the data
#########################################################


#########################################################
#Sub to get data on Informativeness
sub Informative {
	$states = $#SymbolsII + 1;
	@CharStatesSub = ();
	@StateListSub = ();
	@PolyListSub = ();
	$FoundPoly = 0;

	#make starting set
	for($u = 0; $u <= $states; $u++){
		push(@StateListSub, 0)
		}

	#if gap symbol is used, treat appropriately
	if($FoundGap == 1){
		push(@StateListSub, '0');
		}

	#check for different scores for this char
	for($t = 1; $t <= $NCHAR; $t++){
		&LittleInf;
		push(@CharStatesSub, [@StateListSub]);
		push(@PolyListSub, [@TempPoly]);
		}
	}
#end of sub to get data on Informativeness
#########################################################


#########################################################
#sub to get data on Informativeness of single char
sub LittleInf {
	#requires char number $t, on @NexusArraySub
	@TempPoly = ();

	for($u = 0; $u <= $states; $u++){
		$StateListSub[$u] = 0;
		}
	if($FoundGap == 1){
		$StateListSub[$#StateListSub] = 0;
		}

	#check each state
	for($StateCount = 0; $StateCount <= ($states - 1); $StateCount++){
		$CheckSymbol = $SymbolsII[$StateCount];
		for($SpeciesLoop = 0; $SpeciesLoop <= ($NTAX - 1); $SpeciesLoop++){
			$check = $NexusArraySub[$SpeciesLoop][$t];
			if($CheckSymbol eq $check){
				$StateListSub[$StateCount]++
				}
			}
		}

	#check for missing, gap and poly
	for($SpeciesLoop = 0; $SpeciesLoop < $NTAX; $SpeciesLoop++){
		$check = $NexusArraySub[$SpeciesLoop][$t];
		if($Missing eq $check){
			$StateListSub[$#StateListSub]++;
			}
		if($FoundGap == 1){
			if($Gap eq $check){
				$StateListSub[$#StateListSub - 1]++;
				}
			}
		if($check =~ /\)$|\}$/){
			push(@TempPoly, $check);
			$FoundPoly = 1;
			}
		}
	}
#end of sub to get data on Informativeness of single char
#########################################################


#########################################################
#sub to Open new nexus file
sub NewNexusOpening {
	$InterNewNexus = 0;
	open(NEWNEXUS, $OutTypeNex, $PathNexusFile);

	if($Interleaf == 1){
		$InterNewNexus = 1;
		}

# Note PAUP* doesn't seem to able to cope with
# transposed, interleaved matrices!
#	if($Inverted == 1){
#		if($NTAX > 50){
#			$InterNewNexus = 1;
#			}
#		}
	if($Inverted == 0){
		if($NCHAR > 50){
			$InterNewNexus = 1;
			}
		}
	}
#end of sub asking about new matrix
#########################################################


#########################################################
#sub to print interleaved nexus blocks
sub INTERLEAFPRINT {
	@LenghthCheck = @{$NexusArrayINTER[0]};
	$printlength = $#LenghthCheck;

	$InterleafPrintLengthCheck = 1;
	if($printlength < 50){
		$StartChar = 1;
		$PrintNCHAR = $NCHAR;
		@NexusArrayTemp = @NexusArrayINTER;
		&PrintNexus;
		}

	else {
		print $PrintTo "$Explain";
		@tempArray = ();
		$q = 0;
		for($t=1; $t <= $printlength; $t++){
			$q++;
			$CheckNO = $t;
			if((($CheckNO + 50) % 50 == 0) or ($t == $printlength)){
				if($t == $printlength){
					$PrintNCHAR = $t % 50;
					$StartChar = $NCHAR - (($t % 50) - 1)
					}
				else{
					$PrintNCHAR = 50;
					$StartChar = $t - 49
					}

				if($PrintTo eq "NEWNEXUS"){
					print $PrintTo "\n\[$CharTaxRow $StartChar to $t]\n";
					}
				else{
					print $PrintTo "\n<p><br><b>\[$CharTaxRow $StartChar to $t]</b></p>\n";
					}


				for($r = 0; $r < $NTAX; $r++){
					$tempArray[$r][0] = $NexusArrayINTER[$r][0]
					}
				for($r = 0; $r < $NTAX; $r++){
					$tempArray[$r][$q] = $NexusArrayINTER[$r][$t]
					}

				@NexusArrayTemp = @tempArray;
				&PrintNexus;
				$q = 0;
				@tempArray = ();
				}

			else{
				for($r = 0; $r < $NTAX; $r++){
					$tempArray[$r][$q] = $NexusArrayINTER[$r][$t]
					}
				}
			}
		}
	}
#end of sub to print interleaved nexus blocks
#########################################################


#########################################################
#sub to print ordinary nexus blocks
sub PrintNexus {
	#OUTPUT or NEWNEXUS
	#Open html table
	$NBSP = " ";

	if($PrintTo eq "OUTPUT"){
		$NBSP = "&nbsp;";
		print $PrintTo "\n<table class=\"nexus\">";
		}

	#Get longest name for spacing, etc.
	$OriginalNamelength = 0;
	for($i = 0; $i <= $#NexusArrayTemp; $i++){
		@Temp = @{$NexusArrayTemp[$i]};
		foreach $z(@Temp){
			if(length($z) > $OriginalNamelength){
				$OriginalNamelength = length($z)
				}
			}
		}

	if($OriginalNamelength < 5){$OriginalNamelength = 5}
	$Namelength = $OriginalNamelength;
	if($Namelength < 11){$Namelength = 11}

	#Character numbers
	if($PrintTo eq "NEWNEXUS"){
		print $PrintTo "$BrakishO";

		$Shorter = $OriginalNamelength - 1;
		print $PrintTo " "x $Shorter;
		print  $PrintTo "\t";

		$Adjuster = $StartChar - 1;
		for($i = 5; $i <= $PrintNCHAR; $i+=5){
			$j = $Adjuster + $i;
			if($j == 5){
				$Shorter = 0;
				for($Yetti = 0; $Yetti <= 3; $Yetti++){
					$Shorter = $Shorter + $MaxLengths[$Yetti];
					}
				print $PrintTo "$NBSP" x $Shorter;
				$Shorter = $MaxLengths[4] - 1;
				print  $PrintTo "5";
				print $PrintTo "$NBSP" x $Shorter;
				print  $PrintTo "\t";

				}
			else{
				$Shorter = 0;
				for($Yetti = ($j-5); $Yetti <= ($j-2); $Yetti++){
					$Shorter = $Shorter + $MaxLengths[$Yetti];
					}

				print $PrintTo "$NBSP" x $Shorter;
				if(length($j) <= 4){
					$Shorter = $MaxLengths[$j-1] - length($j);
					print $PrintTo "$j";

					}
				else{	$Shorter = $MaxLengths[$j-1] - 1;
					print $PrintTo ".";
					}

				print $PrintTo "$NBSP" x $Shorter;
				print $PrintTo "\t";
				}
			}
		}
	else{
		print $PrintTo "\n<tr>\n\t<td class=\"nexusDat\">&nbsp;</td>\n\t<td class=\"nexusDat\"><div align=\"right\">$BrakishO";
		$Adjuster = $StartChar - 1;
		for($i = 5; $i <= $PrintNCHAR; $i+=5){
			$j = $Adjuster + $i;
			print $PrintTo "$j\&nbsp;</div></td>\n\t<td class=\"nexusDat\"><div align=\"right\">";
			}
		}

	#Character ticks
	if($PrintTo eq "NEWNEXUS"){
		print  $PrintTo "\n";
		print  $PrintTo "$CharTaxCol";
		print $PrintTo " "x(($OriginalNamelength)-5);
		print $PrintTo "\t";

		for($i = 5; $i <= $PrintNCHAR; $i+=5){
			$j = $Adjuster + $i;
			$Shorter = 0;

			if($j == 5){
				for($Yetti = 0; $Yetti <= 3; $Yetti++){
					$Shorter = $Shorter + $MaxLengths[$Yetti];
					}
				print $PrintTo "$NBSP" x $Shorter;
				$Shorter = $MaxLengths[4] - 1;
				print  $PrintTo  ".";
				print $PrintTo "$NBSP" x $Shorter;
				print  $PrintTo  "\t";
				}
			else{
				for($Yetti = ($j-5); $Yetti <= ($j-2); $Yetti++){
					$Shorter = $Shorter + $MaxLengths[$Yetti];
					}
				print $PrintTo "$NBSP" x $Shorter;
				$Shorter = $MaxLengths[$j - 1] - 1;
				print  $PrintTo  ".";
				print $PrintTo " " x $Shorter;
				print  $PrintTo  "\t";
				}
			}
		print  $PrintTo "$BrakishC";
		}
	else{
		print  $PrintTo "\n</div>\n\t</td>\n</tr>\n<tr>";
		print  $PrintTo "<td class=\"nexusDat\"><b>$CharTaxCol</b></td><td class=\"nexusDat\"><div align=\"right\">";
		for($i = 5; $i <= $PrintNCHAR; $i+=5){
			$j = $Adjuster + $i;
			print  $PrintTo  ".&nbsp;</div></td>\n\t<td class=\"nexusDat\"><div align=\"right\">";
			}
		print  $PrintTo "$BrakishC</div></td>\n</tr>\n<tr><td>";
		}

	#Print data nexus with a tab between every 5 characters
	$Dulloff = "";
	for($i = 0; $i <= $#NexusArrayTemp; $i++){
		@Temp = @{$NexusArrayTemp[$i]};
		$x = $#Temp;
		$ListPosition = 0;
		until($ListPosition > $x){
			$z = $Temp[$ListPosition];

			if($PrintTo eq "NEWNEXUS"){
				#print name
				if($ListPosition == 0){
					print $PrintTo "\n$z";
					print $PrintTo "$NBSP"x (($OriginalNamelength)-(length($z)));
					print $PrintTo "\t"
					}

				#print char
				else{
					$Longi = $MaxLengths[$Adjuster + $ListPosition -1];
					$Longi = $Longi - length($z);
					if(($ListPosition > 1) and (($ListPosition - 1)% 5 == 0)){
						print $PrintTo "\t"
						}
					print $PrintTo "$z";
					print $PrintTo "$NBSP" x $Longi;
					}
				}

			else{
				#print name
				if($ListPosition == 0){
					print $PrintTo "$Dulloff</div>\n\t</td>";
					$Dull1 = "";
					$Dulloff = "";
					$Dull1a = "";
					if($DeletionFinalList[$i] == 0){
						$Dull1 = "<p class=\"Grey\">";
						$Dull1a = "<p class=\"GreyNex\">";
						$Dulloff = "</p>";
						}
					$quoteless = $z;
					$quoteless =~ s/^\'//;
					$quoteless =~ s/\'$//;

					print $PrintTo "\n</tr>";
					print $PrintTo "\n<tr>\n\t<td class=\"nexusSp\">$Dull1\<i>$quoteless\&nbsp;</i>$Dulloff";
					print $PrintTo "\n\t</td>";
					print $PrintTo "\n\t<td class=\"nexusDat\"><div align=\"right\">$Dull1a";
					}

				#print char
				else{
					$Dull2a = "";
					$Dull2b = "";
					if(($FinalExSet[$ListPosition + $Adjuster - 1] == 1) and ($DeletionFinalList[$i] == 1)){
						$Dull2a = "<p class=\"GreyNex\">";
						$Dull2b = "</p>";
						}
					$Longi = $MaxLengths[$Adjuster + $ListPosition -1];
					$Longi = $Longi - length($z);
					if(($ListPosition > 1) and (($ListPosition - 1)% 5 == 0)){
						if($DeletionFinalList[$i] == 0){
							$Dull2a = "<p class=\"GreyNex\">";
							$Dull2b = "</p>";
							}
						
						print $PrintTo "&nbsp;$Dulloff\</div>\n\t</td>";
						print $PrintTo "\n\t<td class=\"nexusDat\"><div align=\"right\">$Dull2a";
						}
					if($Dulloff eq '</p>'){
						print $PrintTo "$z";
						print $PrintTo "$NBSP" x $Longi;
						}
					else{
						print $PrintTo "$Dull2a" . "$z";
						print $PrintTo "$NBSP" x $Longi;
						print $PrintTo  "$Dull2b";
						}
					}
				}

			$ListPosition++;
			}
		}
	if($PrintTo eq "OUTPUT"){print $PrintTo "$Dulloff</div>\n\t</td>\n\t</tr></table>\n"}
	}
#end of sub to print normal nexus blocks
#########################################################


#########################################################
#sub to read in polymorphic codings
sub PolymorphReader {
	#print OUTPUT "\n\nI'm afraid that PerlEQ can't yet deal with Multistate codings (i.e. polymorphic taxa)";
	#die "\n\nI'm afraid that PerlEQ can't yet deal with Multistate codings (i.e. polymorphic taxa)\n\n";
	$PolyKill = 0;

	$AreTherePolym = 1;
	$TempPoly = $TempChar;
	$TempChar =  shift(@TempChopper);
	$TempPoly = $TempPoly . $TempChar;

	#read in char, and make into a list
	until($TempChar =~ /\)|\}/){

		if(@TempChopper == 0){
			$StartCol++;
			$LineIn = $UserNEXUSstore[$StartCol];
			unless($RespectCase == 1){
				$LineIn =~ tr/a-z/A-Z/;
				}
			if($LineIn =~ /\n$/){
				$Breaker = 1;
				}
			chomp($LineIn);
			@TempChopper = split(//, $LineIn);
			}

		$TempChar =  shift(@TempChopper);
		$TempPoly = $TempPoly . $TempChar;
		}

	$TempChar = "";
	$TempPoly =~ s/,//;
	@TempPolyList = split(//, $TempPoly);

	#check each code for gap, matchchar, unknown, etc
	if($DoIcheck == 1){
		for($th = 1; $th <= ($#TempPolyList - 1); $th++){
			$TempChar = $TempPolyList[$th];
			unless($TaxaCountDown == $NTAX){
				unless($Matchchar eq ""){	
					if("$TempChar" eq "$Matchchar"){
						if($Inverted == 1){
							if(@SpeciesData == 0){
								$TempChar = $TempList[0];
								}
							else{
								$TempChar = $SpeciesData[0];
								}
							}
						else{
							$TempChar = $NexusII[0][($NCHAR - $CharCountDown)];
							};
						}
					}	
				}

			&CharacterCheck;
			if($GoodChar == 0){
				push(@UnrecList, $UnrecString);
				$TempChar = "";
				}
			$TempPolyList[$th] = $TempChar;

			}
		}

	#reconstruct as text, in order of symbols list
	@TempPolyList = sort{&PolySorter}(@TempPolyList);
	$TempChar = "";
	if(@TempPolyList == 3){
		$TempChar = $TempPolyList[1];
		}
	else{
		foreach $GH (@TempPolyList){
			unless($TempChar =~ quotemeta($GH)){
				$TempChar .= $GH;
				}

			}
		}

	#convert bracket types, if specified
	if($MStaxa eq "polymorph"){
		$TempChar =~ s/\{/\(/;
		$TempChar =~ s/\}/\)/;
		}
	if($MStaxa eq "uncertain"){
			$TempChar =~ s/\(/\{/;
			$TempChar =~ s/\)/\}/;
		}

	#Temp patch to kill polymorphic scores
	if($TempChar =~ /\)/){
		$TempChar = $Missing;
		$PolyKill = 1;
		}

	}
#end of sub to read in polymorphic codings
#########################################################


#########################################################
#sub to get max char lengths for printing
sub MaxCharLengths {
	#reset lenghts
	@MaxLengths = ();
	for($i = 1; $i <= $NCHAR; $i++){
		push(@MaxLengths, '1')
		}

	#check each char
	for($i = 1; $i <= $NCHAR; $i++){
		for($s = 0; $s <= $#LengthNexusArray; $s++){
			if($MaxLengths[$i-1] < length($LengthNexusArray[$s][$i])){
				$MaxLengths[$i-1] = length($LengthNexusArray[$s][$i]);
				}
			}
		}
	}
#end of sub to get max char lengths for printing
#########################################################


#########################################################
#tiny sub to sort polymorphic codes into order
sub PolySorter {
	#1 = $b,$a	-1 = $a,$b
	if($a =~ /\(|\{/){
		return(-1)
		}
	elsif($a =~ /\)|\}/){
		return(1)
		}
	elsif($a eq $Gap){
		if($b =~ /\)|\}/){
			return(-1)
			}
		else{
			return(1)
			}
		}

	else{
		foreach $Bit (@SymbolsII){
			if($a eq $Bit){
				if($b =~ /\(|\{/){
					return(1)
					}
				else{
					return(-1)
					}
				last
				}
			elsif($b eq $Bit){
				if($a =~ /\(|\{/){
					return(-1)
					}
				else{
					return(1)
					}
				}
			}
		return(0)
		}
	}
#end of tiny sub to sort polymorphic codes into order
#########################################################


#########################################################
#sub to read ordinary Nexus files
sub ReadNexus {
	$TaxaCountDown = $NTAX;
	$StartCol++;
	$FromMatrixTest = 0;
	$DoIcheck = 1;

	#Big LOOP to get all the nexus data...
	while($TaxaCountDown >= 1){
		@SpeciesData = ();
		$CharCountDown = $NCHAR;

		#get species names, is first
		if($MatrixLabels eq "left"){
			$SpeciesName = $UserNEXUSstore[$StartCol];
			chomp($SpeciesName);

			if($Inverted == 1){
				$FromMatrixChars[$NTAX - $TaxaCountDown] = $SpeciesName;
				$FromMatrixCharTest = 1;
				}
			if($SpeciesName ne "_"){
				$FromMatrixSpecies[$NTAX - $TaxaCountDown] = $SpeciesName;
				$FromMatrixTest = 1;
				}
			$StartCol++;
			}


		#get following data
		until($CharCountDown < 1){
			$LineIn = $UserNEXUSstore[$StartCol];
			unless($RespectCase == 1){
				$LineIn =~ tr/a-z/A-Z/;
				}
			chomp($LineIn);
			@TempList = ();
			@TempChopper = split(//, $LineIn);
			while(@TempChopper > 0){
				$TempChar = shift(@TempChopper);

				if($TempChar =~ /\(|\{/){
					&PolymorphReader;
					unless(($TaxaCountDown == $NTAX) or ($Matchchar eq "")){
						if("$TempChar" eq "$Matchchar"){
							&MatchingChars;
							}
						}
					&CharacterCheck;
					unless(($TaxaCountDown == $NTAX) or ($Matchchar eq "")){
						if("$TempChar" eq "$Matchchar"){
							&MatchingChars;
							}
						}
					$CharCountDown--;
					push(@TempList, $TempChar);
					next;
					}

				unless(($TaxaCountDown == $NTAX) or ($Matchchar eq "")){
					if("$TempChar" eq "$Matchchar"){
						&MatchingChars;
						}
					}
				&CharacterCheck;
				unless(($TaxaCountDown == $NTAX) or ($Matchchar eq "")){
					if("$TempChar" eq "$Matchchar"){
						&MatchingChars;
						}
					}

				if($GoodChar == 0){
					push(@UnrecList, $UnrecString);
					$TempChar = $Missing;
					}
				$CharCountDown--;
				push(@TempList, $TempChar)
				}
			push(@SpeciesData, @TempList);
			$StartCol++;
			}

		#get species names, if last
		if($MatrixLabels eq "right"){
			$SpeciesName = $UserNEXUSstore[$StartCol];
			chomp($SpeciesName);

			if($Inverted == 1){
				$FromMatrixChars[$NTAX - $TaxaCountDown] = $SpeciesName;
				$FromMatrixCharTest = 1;
				}
			if($SpeciesName ne "_"){
				$FromMatrixSpecies[$NTAX - $TaxaCountDown] = $SpeciesName;
				$FromMatrixTest = 1;
				}
			$StartCol++;
			}

		push(@NexusII, [@SpeciesData]);
		$TaxaCountDown--;
		}
	@NexusArray = @NexusII;
	}
#end of sub to read ordinary Nexus files
#########################################################


#########################################################
#sub to read interleafed Nexus files
sub ReadInterLeafNexus {
	$DoIcheck = 1;
	$TaxaCountDown = $NTAX;
	$StartCol++;
	#Big LOOP to get all the nexus data...
	while($TaxaCountDown >= 1){
		@SpeciesData = ();
		$CharCountDown = $NCHAR;

		#get species names, if first
		if($MatrixLabels eq "left"){
			$SpeciesName = $UserNEXUSstore[$StartCol];
			chomp($SpeciesName);
			if($Inverted == 1){
				$FromMatrixChars[$NTAX - $TaxaCountDown] = $SpeciesName;
				$FromMatrixCharTest = 1;
				}
			if($SpeciesName ne "_"){
				$FromMatrixSpecies[$NTAX - $TaxaCountDown] = $SpeciesName;
				$FromMatrixTest = 1;
				}
			$StartCol++;
			}

		#Read in first line of data
		$LineIn = $UserNEXUSstore[$StartCol];
		unless($RespectCase == 1){
			$LineIn =~ tr/a-z/A-Z/;
			}
		$Breaker = 0;
		until($Breaker == 1){
			if($LineIn =~ /\n$/){
				$Breaker = 1;
				$ColumnStore = ($StartCol + 1);
				if($MatrixLabels eq "right"){
					chomp($LineIn);
					if($Inverted == 1){
						$FromMatrixChars[$NTAX - $TaxaCountDown] = $SpeciesName;
						$FromMatrixCharTest = 1;
						}
					if($LineIn ne "_"){
						$FromMatrixSpecies[$NTAX - $TaxaCountDown] = $LineIn;
						$FromMatrixTest = 1;
						}
					$StartCol++;
					$LineIn = $UserNEXUSstore[$StartCol];
					unless($RespectCase == 1){
						$LineIn =~ tr/a-z/A-Z/;
						}
					last
					}
				}
			chomp($LineIn);
			$LineIn =~ s/\s//;

			@TempList = ();
			@TempChopper = split(//, $LineIn);
			while(@TempChopper > 0){
				$TempChar =  shift(@TempChopper);

				if($TempChar =~ /\(|\{/){
					&PolymorphReader;
					unless(($TaxaCountDown == $NTAX) or ($Matchchar eq "")){
						if("$TempChar" eq "$Matchchar"){
							&MatchingChars;
							}
						}
					&CharacterCheck;
					unless(($TaxaCountDown == $NTAX) or ($Matchchar eq "")){
						if("$TempChar" eq "$Matchchar"){
							&MatchingChars;
							}
						}
					$CharCountDown--;
					push(@TempList, $TempChar);
					next;
					}

				unless(($TaxaCountDown == $NTAX) or ($Matchchar eq "")){
					if("$TempChar" eq "$Matchchar"){
						&MatchingChars;
						}
					}
				&CharacterCheck;
				unless(($TaxaCountDown == $NTAX) or ($Matchchar eq "")){
					if("$TempChar" eq "$Matchchar"){
						&MatchingChars;
						}
					}
				if($GoodChar == 0){
					push(@UnrecList, $UnrecString);
					$TempChar = $Missing
					}
				$CharCountDown--;
				push(@TempList, $TempChar)
				}
			push(@SpeciesData, @TempList);
			$StartCol++;
			$LineIn = $UserNEXUSstore[$StartCol];
			unless($RespectCase == 1){
				$LineIn =~ tr/a-z/A-Z/;
				}
			}
		@TempList = ();

		#Look for more rows of data
		#while loop to get rest of data for that species
		while($CharCountDown >= 1){
			#find next
			&FindNextBlock;

			#jump species names, if first
			if($MatrixLabels eq "left"){
				$StartCol++;
				}

			#Read in next line of data
			$LineIn = $UserNEXUSstore[$StartCol];
			unless($RespectCase == 1){
				$LineIn =~ tr/a-z/A-Z/;
				}
			$Breaker = 0;
			until($Breaker == 1){
				if($LineIn =~ /\n$/){
					$Breaker = 1;
					if($MatrixLabels eq "right"){
						$StartCol++;
						$LineIn = $UserNEXUSstore[$StartCol];
						unless($RespectCase == 1){
							$LineIn =~ tr/a-z/A-Z/;
							}
						last
						}
					};
				chomp($LineIn);
				$LineIn =~ s/\;//;

				@TempList = ();
				@TempChopper = split(//, $LineIn);
				while(@TempChopper > 0){
					$TempChar =  shift(@TempChopper);

					if($TempChar =~ /\(|\{/){
						&PolymorphReader;
						unless(($TaxaCountDown == $NTAX) or ($Matchchar eq "")){
							if("$TempChar" eq "$Matchchar"){
								&MatchingChars;
								}
							}
						&CharacterCheck;
						unless(($TaxaCountDown == $NTAX) or ($Matchchar eq "")){
							if("$TempChar" eq "$Matchchar"){
								&MatchingChars;
								}
							}
						$CharCountDown--;
						push(@TempList, $TempChar);
						next;
						}

					unless(($TaxaCountDown == $NTAX) or ($Matchchar eq "")){
						if("$TempChar" eq "$Matchchar"){
							&MatchingChars;
							}
						}
					&CharacterCheck;
					unless(($TaxaCountDown == $NTAX) or ($Matchchar eq "")){
						if("$TempChar" eq "$Matchchar"){
							&MatchingChars;
							}
						}

					if($GoodChar == 0){
						push(@UnrecList, $UnrecString);
						$TempChar = $Missing
						}
					$CharCountDown--;
					push(@TempList, $TempChar)
					}
				push(@SpeciesData, @TempList);
				$StartCol++;
				$LineIn = $UserNEXUSstore[$StartCol];
				unless($RespectCase == 1){
					$LineIn =~ tr/a-z/A-Z/;
					}
				}
			@TempList = ();
			if($StartCol > $#UserNEXUSstore){last};
			}# end of while loop to get rest of data for that species

		$TaxaCountDown--;
		push(@NexusII, [@SpeciesData]);
		$StartCol = $ColumnStore;
		if($StartCol > $#UserNEXUSstore){last};
		}
	@NexusArray = @NexusII;
	}
#end of sub to read interleafed Nexus files
#########################################################


#########################################################
#sub to replace matchchar symbol with real code
sub MatchingChars {
	if($Inverted == 1){
		if(@SpeciesData == 0){
			$TempChar = $TempList[0];
			}
		else{
			$TempChar = $SpeciesData[0];
			}
		}
	else{
		$TempChar = $NexusII[0][($NCHAR - $CharCountDown)];
		};
	}
#end of sub to replace matchchar symbol with real code
#########################################################


#########################################################
#little sub to find next block
sub FindNextBlock {
	$StartCol--;
	$LineIn = $UserNEXUSstore[$StartCol];
	$CycleRound = $NTAX;
	while($CycleRound > 0){
		if($LineIn =~ /\n$/){
			$CycleRound--
			}
		$StartCol++;
		$LineIn = $UserNEXUSstore[$StartCol];
		}
	}
#end of little sub to find next block
#########################################################


#########################################################
#Sub to check characters to see if they match the Symbols list or Missing
sub CharacterCheck {
	$GoodChar = 0;

	foreach $pp (@SymbolsListII){
		if($TempChar eq $pp){$GoodChar = 1}
		}
	unless(@From == 0){
		for($pp = 0; $pp <= $#From; $pp++){
			if($TempChar eq $From[$pp]){
				$TempChar = $To[$pp];
				$GoodChar = 1;
				};
			}
		}


	if($GoodChar == 0){
		$Doog = 1 + $NCHAR - $CharCountDown;
		$UnrecString = "<i\>$SpeciesName\<\/i\> character $Doog\: symbol \'$TempChar\' unrecognised, changed to \'$Missing\'";
		}

	}
#End of sub to check characters
#########################################################


#########################################################
#Sub to create various output files
sub OpenOut {
	while(1){
		print "\nSpecify or create the name for the output files\n";
		print "The program will save a data \(.htm\) and a nexus \(.nex\) file\n";
		print "Filenames of more than 25 characters will be shortened:\n   ";

		$OutFile = <stdin>;
		chomp($OutFile);
		@FileDetails = fileparse($OutFile, @SuffixList);
		$OtherPath = $FileDetails[1];
		if($OtherPath eq $LocalPath){
			$OtherPath = $MyPath;
			}
		$OutFile = $FileDetails[0];


		$lengthtest=0;

		# make sure output file is of the right length
		$OldOutFile = $OutFile;

		if($OutFile eq "q"){die "\nProgram quit"}
		elsif(length($OutFile) > 25){
			until(length($OutFile) <= 25){chop $OutFile}
			$lengthtest=1;
			last
			}

		else{last}
		}#end of while

	$DataFile = $OutFile . ".htm";
	$PathDataFile = $OtherPath . $DataFile;
	$NexusFile = $OutFile . ".nex";
	$PathNexusFile = $OtherPath . $NexusFile;

	if ($lengthtest == 1){
		print "\n$OldOutFile shortened...\nOutput file saved as $OutFile\n"
		}

	#check if out file exists
	$loopend = 0;
	$OutTypeDat = ">>";
	$OutTypeNex = ">>";

	while(1){

		if (-e $PathDataFile){
			print "\nThe data file $DataFile already exists.  Append (a) or Replace (r)\n   ";
			$Exists = <stdin>;
			chomp($Exists);
			$Exists =~ tr/A-Z/a-z/;
			if($Exists eq "a"){
				$OutTypeDat = ">>";
				last
				}
			elsif($Exists eq "r"){
				$OutTypeDat = ">";
				last
				}
			else {next}
			}
		last}

	while(1){
		if ((-e $PathNexusFile) != 0){
			print "\nThe nexus file $NexusFile already exists.  Append (a) or Replace (r)\n   ";
			$Exists = <stdin>;
			chomp($Exists);
			$Exists =~ tr/A-Z/a-z/;
			if($Exists eq "a"){
				$OutTypeNex = ">>";
				last
				}
			elsif($Exists eq "r"){
				$OutTypeNex = ">";
				last
				}
			else {next}
			}
		last}

	}
#end of sub to open output file
#########################################################


#########################################################
#sub to find and record PAUP blocks
sub PAUPBLOCK {
	@PAUPblockTotal = ();
	@ExcludeTaxaNEXUS = @UserNEXUSstore;

	@DeletionFinalList = ();
	for($JJK = 0; $JJK < $NTAX; $JJK++){
		$Adder = 1;
		push(@DeletionFinalList, $Adder);
		}

	#Read in any PAUP blocks
	$BeginTest = 0;
	for($fk = 0; $fk <= $#ExcludeTaxaNEXUS; $fk++){
		$ExcludedTaxaLine = $ExcludeTaxaNEXUS[$fk];
		chomp($ExcludedTaxaLine);
		if($ExcludedTaxaLine =~ /^begin$/i){
			$BeginTest = 1;
			next;
			}
		if(($BeginTest == 1) and ($ExcludedTaxaLine =~ /^paup$/i)){
			&ReadInPAUPBlock;
			$BeginTest = 0;
			next;
			}
		else{$BeginTest = 0}
		}#end of search through whole Nexus file for PAUP blocks
	}
#sub to find and record PAUP blocks
#########################################################


#########################################################
#sub to exclude taxa
sub ExcludeTaxa {
	#Now check for Delete commands
	for($LookingAtPAUP = 0; $LookingAtPAUP <= $#PAUPblockTotal; $LookingAtPAUP++){
		$TestPAUP = $PAUPblockTotal[$LookingAtPAUP];
		if($TestPAUP =~ /^mstaxa$/i){
			$LookingAtPAUP += 2;
			$TestPAUP = $PAUPblockTotal[$LookingAtPAUP];
			$TestPAUP =~ tr/A-Z/a-z/;
			$MStaxa = $TestPAUP;
			$FoundMStaxa = 1;
			}
		if($TestPAUP =~ /^delete$/i){&DeleteProcedure};
		}
	}
#end of sub excluded taxa
#########################################################


#########################################################
#sub to read in PAUP blocks
sub ReadInPAUPBlock {
	$fk++;
	$TempTestLine = $ExcludeTaxaNEXUS[$fk];
	chomp($TempTestLine);

	until($TempTestLine =~ /^end$/i){
		push(@PAUPblockTotal, $TempTestLine);
		$fk++;
		$TempTestLine = $ExcludeTaxaNEXUS[$fk];
		chomp($TempTestLine);
		}#end of until

	$AreThereTaxSets = 1;
	}
#end of sub to read in PAUP blocks
#########################################################


#########################################################
#sub to exclude characters in PAUP block
sub ExcludeFinder {
	for($LookingAtPAUP = 0; $LookingAtPAUP <= $#PAUPblockTotal; $LookingAtPAUP++){
		$TestPAUP = $PAUPblockTotal[$LookingAtPAUP];
		if($TestPAUP =~ /^exclude$/i){
			&ExcludeProcedure
			};
		}

	}
#end of sub to exclude characters in PAUP block
#########################################################


#########################################################
#sub to actually exclude the PAUP block characters
sub ExcludeProcedure {
	$SomeExcluded = 0;

	#find species names or labels
	until($TestPAUP =~ /^;$/){
		$LookingAtPAUP++;

		$TestPAUP = $PAUPblockTotal[$LookingAtPAUP];
		$NextN = $PAUPblockTotal[1 + $LookingAtPAUP];
		chomp($TestPAUP);
		chomp($NextN);
		if($TestPAUP =~ /^\/$/){next};
		if($TestPAUP =~ /^only$/){next};

		@LookingSet = @{$TotalCharLabels{"$TestPAUP"}};

		#if a charset
		if($LookingSet[0] eq "charset"){

			shift(@LookingSet);
			$CharCountered = 0;
			foreach $rw (@LookingSet){
				if($rw == 1){
					$FinalExSet[$CharCountered] = 1;
					$SomeExcluded = 1;
					}
				$CharCountered++;
				}
			next;
			}

		#check for a range
		elsif($NextN eq "-"){
			$LookingAtPAUP += 2;
			$SndTestPAUP = $PAUPblockTotal[$LookingAtPAUP];
			@SecondSet = @{$TaxLabelTotalArray{"$SndTestPAUP"}};

			$IncreNot = 0;
			$Incrementer = $PAUPblockTotal[$LookingAtPAUP + 1];
			if($Incrementer =~ /^\\$/){
				$Incrementer = $PAUPblockTotal[$LookingAtPAUP + 2];
				if($Incrementer =~ /^\d$/){
					$IncreNot = 1;
					}
				else{
					$Incrementer = 1;
					$IncreNot = 0;
					}
				}
			else{$Incrementer = 1;}


			#if all
			if($SndTestPAUP eq "\."){
				for($DFG = $LookingSet[1]; ($DFG - 1) <= $#FinalExSet; $DFG += $Incrementer){

					$FinalExSet[$DFG - 1] = 1;
					$SomeExcluded = 1;
					}
				}#end of if

			else{
				for($DFG = $LookingSet[1]; $DFG <= $SecondSet[1]; $DFG += $Incrementer){
					$FinalExSet[$DFG - 1] = 1;
					$SomeExcluded = 1;
					}
				}#end of else
			if($IncreNot == 1){$LookingAtPAUP += 2}
			next;
			}#end of small if to check for a range

		#if a singleton
		elsif(($LookingSet[0] eq "char") and ($NextN ne "-")){
			$FinalExSet[$LookingSet[1] - 1] = 1;
			$SomeExcluded = 1;
			}
		@LookingSet = ();
		@SecondSet = ();
		}

	#check eliminate statements
	if($SomeExcluded == 1){
		$AllIncTest = 0;
		$SomeEliminated = 0;
		&EliminateStatement; #gives @EliminateList, @EliminateDiff & mod @FinalExSet
		}

	if(($SomeExcluded == 1) or ($SomeEliminated == 1)){
		for($GT = 0; $GT <= $#NexusArray; $GT++){
			for($ED = 1; $ED <= $NCHAR; $ED++){
				if($FinalExSet[$ED-1] == 1){
					$NexusArray[$GT][$ED] = $Missing
					}
				}
			}
		}
	}
#end of sub to actually exclude the PAUP block characters
#########################################################


#########################################################
#sub to actually delete taxa
sub DeleteProcedure {
	$SomeDeleted = 0;

	#find species names or labels
	until($TestPAUP =~ /^;$|^\/$/){
		$LookingAtPAUP++;
		$TestPAUP = $PAUPblockTotal[$LookingAtPAUP];
		$NextN = $PAUPblockTotal[1 + $LookingAtPAUP];
		chomp($TestPAUP);
		chomp($NextN);

		@LookingSet = @{$TaxLabelTotalArray{"$TestPAUP"}};

		#if a taxset
		if($LookingSet[0] eq "taxset"){
			shift(@LookingSet);
			$TaxonCountered = 0;
			foreach $rw (@LookingSet){
				if($rw == 1){
					$DeletionFinalList[$TaxonCountered] = 0;
					$SomeDeleted = 1;
					}
				$TaxonCountered++;
				}
			next;
			}

		#check for a range
		elsif($NextN eq "-"){
			$LookingAtPAUP += 2;
			$SndTestPAUP = $PAUPblockTotal[$LookingAtPAUP];
			@SecondSet = @{$TaxLabelTotalArray{"$SndTestPAUP"}};

			#if all
			if($SndTestPAUP eq "\."){
				for($DFG = $LookingSet[1]; $DFG <= $#DeletionFinalList; $DFG++){
					$DeletionFinalList[$DFG - 1] = 0;
					$SomeDeleted = 1;
					}
				}#end of if

			else{
				for($DFG = $LookingSet[1]; $DFG <= $SecondSet[1]; $DFG++){
					$DeletionFinalList[$DFG - 1] = 0;
					$SomeDeleted = 1;
					}
				}#end of else
			next;
			}#end of small if to check for a range

		#if a singleton
		elsif(($LookingSet[0] eq "taxon") and ($NextN ne "-")){
			$DeletionFinalList[$LookingSet[1] - 1] = 0;
			$SomeDeleted = 1;
			}
		@LookingSet = ();
		@SecondSet = ();
		}

	if($SomeDeleted == 1){
		for($DeleteThem = 0; $DeleteThem <= $#NexusArray; $DeleteThem++){
			if($DeletionFinalList[$DeleteThem] == 0){
				for($MakeQuest = 1; $MakeQuest <= $NCHAR; $MakeQuest++){
					$NexusArray[$DeleteThem][$MakeQuest] = $Missing;
					}
				}
			}
		}
	}
#END OF sub to actually delete taxa
#########################################################


#########################################################
# Sub to read in Character sets types
sub FindCharacterSets {
	@CharsetNEXUS = @UserNEXUSstore;

	#make default set 'all'
	@Slider = ('charset');
	for($d = 1; $d <= $NCHAR; $d++){
		push(@Slider, 1);
		}
	$TotalCharLabels{"all"} = [@Slider];

	#make default set 'uninf', 'constant', 'gapped' and 'missambig'
	@Slider = ('charset');
	for($d = 1; $d <= $NCHAR; $d++){
		push(@Slider, 0);
			}
	$TotalCharLabels{"uninf"} = [@Slider];
	$TotalCharLabels{"constant"} = [@Slider];
	$TotalCharLabels{"gapped"} = [@Slider];
	$TotalCharLabels{"missambig"} = [@Slider];

	for($fk = 0; $fk <= $#CharsetNEXUS; $fk++){
		$BT = $CharsetNEXUS[$fk];
		chomp($BT);
		&LookAroundChars;
		}#end of search through whole Nexus file
	}
#end sub to find Character Sets
#########################################################


#########################################################
#Sub to look through list for charsets
sub LookAroundChars {
	$StandardCheck = 0;
	if($BT =~ /^charset$/i){
		$fk++;
		$CharsetName = $CharsetNEXUS[$fk];
		chomp($CharsetName);

		if($CharsetName eq "\*"){
			$fk++;
			$CharsetName = $CharsetNEXUS[$fk];
			chomp($CharsetName);
			}

		@ThisCharSet = ('charset');
		$fk++;
		$CharsetType = $CharsetNEXUS[$fk];
		$CharsetType =~ tr/A-Z/a-z/;

		####Deal with vector
		if($CharsetType =~ /^\(vector\)$/i){
			$StandardCheck = 0;
			$fk += 2;
			$LookingLine = $CharsetNEXUS[$fk];
			chomp($LookingLine);
			until($LookingLine =~ /;/){
				@TempLooking = ();
				until(length($LookingLine) == 0){
					$LaLaLa = chop($LookingLine);
					unshift(@TempLooking, $LaLaLa);
					}
				push(@ThisCharSet, @TempLooking);
				$fk++;
				$LookingLine = $CharsetNEXUS[$fk];
				chomp($LookingLine);
				}

			$CharSetFinalList{"$CharsetName"} = [@ThisCharSet];
			$TotalCharLabels{"$CharsetName"} = [@ThisCharSet];
			$ThisCharSet[0] = $CharsetName;
			push(@CharsetOnlyList, [@ThisCharSet]);
			$AreThereCharSets = 1;
			@ThisCharSet = ();
			$fk--;
			}

		####Deal with standard
		elsif($CharsetType =~ /^\(standard\)$/i){
			$fk += 2;
			$StandardCheck = 1
			}
		else{
			$fk++;
			$StandardCheck = 1;
			}

		if($StandardCheck == 1){
			for($hG = 1; $hG <= $NCHAR; $hG++){
				push(@ThisCharSet, "0");
				}

			&LookingCharLiner;
			$fk++;
			&SecondCharLooker;

			until($LookingLine =~ /\;/){

				if($LookingLine =~ /^\/$/){next};
				if($LookingLine =~ /^only$/){next};

				#if the first thing is a set
				if($TempLiner[0] eq "charset"){
					shift(@TempLiner);
					$TempNumber = 1;
					foreach $c (@TempLiner){
						if($c == 1){
							$ThisCharSet[$TempNumber] = 1;
							}
						$TempNumber++;
						}
					$fk++;
					&LookingCharLiner;
					$fk++;
					&SecondCharLooker;
					next;
					}

				#else get the char number
				if($TempLiner[0] eq "char"){
					$TempNumber = $TempLiner[1];
					}

				#if there\'s a range
				if($SecondLook eq "-"){
					$fk++;
					&SecondCharLooker;

					$IncreNot = 0;
					$Incrementer = $CharsetNEXUS[$fk + 1];
					if($Incrementer =~ /^\\$/){
						$Incrementer = $CharsetNEXUS[$fk + 2];
						if($Incrementer =~ /^\d$/){
							$IncreNot = 1;
							}
						else{
							$Incrementer = 1;
							$IncreNot = 0;
							}
						}
					else{$Incrementer = 1;}


					#if all
					if($SecondLook =~ /\./){
						for($YH = $TempNumber; $YH <= $NCHAR; $YH += $Incrementer){
							$ThisCharSet[$YH] = 1;
							}
						$fk++;
						&LookingCharLiner;
						$fk++;
						&SecondCharLooker;
						}

					#if a limited range
					else{
						$TempNumberTwo = $SecondTempLiner[1];
						for($YH = $TempNumber; $YH <= $TempNumberTwo; $YH += $Incrementer){
							$ThisCharSet[$YH] = 1;
							}
						}
					if($IncreNot == 1){$fk += 2}

					$fk++;
					&LookingCharLiner;
					$fk++;
					&SecondCharLooker;
					}

				#It second look is a semicolon
				elsif($SecondLook eq "\;"){

					$TempNumber = $TempLiner[1];
					$ThisCharSet[$TempNumber] = 1;
					$LookingLine = $SecondLook;
					}

				#Only a singleton
				else{
					$TempNumber = $TempLiner[1];
					$ThisCharSet[$TempNumber] = 1;
					@TempLiner = @SecondTempLiner;
					$LookingLine = $SecondLook;
					$fk++;
					&SecondCharLooker;
					}#end of else
				}#end of search until ';'

			$CharSetFinalList{"$CharsetName"} = [@ThisCharSet];
			$TotalCharLabels{"$CharsetName"} = [@ThisCharSet];
			$ThisCharSet[0] = $CharsetName;
			push(@CharsetOnlyList, [@ThisCharSet]);
			$StandardCheck = 0;
			$AreThereCharSets = 1;
			$fk--;
			@ThisCharSet = ();
			}#end of input of standard format

		}#end of found charset
	}#end of search through this line of the Nexus file
#end of sub to look through list for charsets
#########################################################


#########################################################
#little sub to get $LookingLine
sub LookingCharLiner {
	@TempLiner = ();
	$LookingLine = $CharsetNEXUS[$fk];
	chomp($LookingLine);

	#check for default sets...
	if($LookingLine =~ /^all$/i){
		$LookingLine = "all";
		}
	if($LookingLine =~ /^uninf$/i){
		$LookingLine = "uninf";
		}
	if($LookingLine =~ /^constant$/i){
			$LookingLine = "constant";
		}
	if($LookingLine =~ /^gapped$/i){
			$LookingLine = "gapped";
		}
	if($LookingLine =~ /^missambig$/i){
				$LookingLine = "missambig";
		}

	#get charset or char
	@TempLiner = @{$TotalCharLabels{"$LookingLine"}};
	}
#end of little sub to get $LookingCharLine
#########################################################


#########################################################
#little sub to get $SecondLook
sub SecondCharLooker {
	@SecondTempLiner =();
	$SecondLook  = $CharsetNEXUS[$fk];
	chomp($SecondLook);
	@SecondTempLiner = @{$TotalCharLabels{"$SecondLook"}};
	}
#end of little sub to get $SecondCharLook
#########################################################


#########################################################
# Sub to read in Taxa sets types
sub FindTaxonSets {
	@TaxsetNEXUS = @UserNEXUSstore;
	@FinalNameList = (); #names on matrix

	#get names on matrix
	for($kj = 0; $kj <= $#NexusArray; $kj++){
		$OTUname = $NexusArray[$kj][0];
		push(@FinalNameList, $OTUname);
		}

	for($fk = 0; $fk <= $#TaxsetNEXUS; $fk++){
		$BT = $TaxsetNEXUS[$fk];
		chomp($BT);
		&LookAroundTaxa;
		}#end of search through whole Nexus file
	}
#end sub FindTaxonSets
#########################################################


#########################################################
#Sub to look through list for taxsets
sub LookAroundTaxa {
	$StandardCheck = 0;
	if($BT =~ /taxset/i){
		$fk++;
		$TaxsetName = $TaxsetNEXUS[$fk];
		chomp($TaxsetName);

		if($TaxsetName eq "\*"){
			$fk++;
			$TaxsetName = $TaxsetNEXUS[$fk];
			chomp($TaxsetName);
			}

		@ThisTaxSet = ('taxset');
		$SetNameTemp = $TaxsetName;
		$fk++;
		$TaxsetName = $TaxsetNEXUS[$fk];

		####Deal with vector
		if($TaxsetName =~ /^\(vector\)$/i){
			$StandardCheck = 0;
			$fk += 2;
			$LookingLine = $TaxsetNEXUS[$fk];
			chomp($LookingLine);
			until($LookingLine =~ /;/){
				@TempLooking = ();
				until(length($LookingLine) == 0){
					$LaLaLa = chop($LookingLine);
					unshift(@TempLooking, $LaLaLa);
					}
				push(@ThisTaxSet, @TempLooking);
				$fk++;
				$LookingLine = $TaxsetNEXUS[$fk];
				chomp($LookingLine);
				}

			$TaxLabelTotalArray{"$SetNameTemp"} = [@ThisTaxSet];
			$ThisTaxSet[0] = $SetNameTemp;
			push(@TaxSetFinalList, [@ThisTaxSet]);
			@ThisTaxSet = ();
			$fk--;
			}

		####Deal with standard
		elsif($TaxsetName =~ /^\(standard\)$/i){
			$fk += 2;
			$StandardCheck = 1
			}
		else{
			$fk++;
			$StandardCheck = 1;
			}

		@LookingSet =();
		@SecondSet =();

		if($StandardCheck == 1){
			for($hG = 1; $hG <= $NTAX; $hG++){
				push(@ThisTaxSet, "0");
				}
			&LookingLiner;
			$fk++;
			&SecondLooker;

			until($LookingLine =~ /\;/){
				#if there\'s a range
				if($SecondLook eq "-"){
					$fk++;
					&SecondLooker;
					#if all
					if($SecondLook =~ /\./){
						for($YH = $LookingSet[1]; $YH <= $NTAX; $YH++){
							$ThisTaxSet[$YH] = 1;
							}
						$fk++;
						&LookingLiner;
						$fk++;
						&SecondLooker;
						}

					#if a limited range
					else{
						for($YH = $LookingSet[1]; $YH <= $SecondSet[1]; $YH++){
							$ThisTaxSet[$YH] = 1;
							}
						}
					}

				#It a semicolon
				elsif($SecondLook eq "\;"){
					if($LookingSet[0] eq "taxset"){
						shift(@LookingSet);
						$Taxcountry = 1;
						foreach $hj (@LookingSet){
							if($hj == 1){
								$ThisTaxSet[$Taxcountry] = 1;
								}
							$Taxcountry++;
							}
						}

					if($LookingSet[0] eq "taxon"){
						$ThisTaxSet[$LookingSet[1]] = 1;
						}
					$LookingLine = $SecondLook;
					}

				#Only a name
				else{
					if($LookingSet[0] eq "taxset"){
						shift(@LookingSet);
						$Taxcountry = 1;
						foreach $hj (@LookingSet){
							if($hj == 1){
								$ThisTaxSet[$Taxcountry] = 1;
								}
							$Taxcountry++;
							}
						}

					if($LookingSet[0] eq "taxon"){
						$ThisTaxSet[$LookingSet[1]] = 1;
						}

					$LookingLine = $SecondLook;
					@LookingSet = @SecondSet;
					$fk++;
					&SecondLooker;
					}#end of else
				}#end of search until ';'

			$TaxLabelTotalArray{"$SetNameTemp"} = [@ThisTaxSet];
			$ThisTaxSet[0] = $SetNameTemp;
			push(@TaxSetFinalList, [@ThisTaxSet]);
			@ThisTaxSet = ();
			$StandardCheck = 0;
			$AreThereTaxSets = 1;
			$fk--;
			}#end of input of standard format
		}#end of found taxset
	}#end of search through this line of the Nexus file

#end of sub to look through list for taxsets
#########################################################


#########################################################
#little sub to get $LookingLine
sub LookingLiner {
	$LookingLine = $TaxsetNEXUS[$fk];
	chomp($LookingLine);
	@LookingSet = @{$TaxLabelTotalArray{"$LookingLine"}};
	}
#end of little sub to get $LookingLine
#########################################################


#########################################################
#little sub to get $SecondLook
sub SecondLooker {
	$SecondLook  = $TaxsetNEXUS[$fk];
	chomp($SecondLook);
	@SecondSet = @{$TaxLabelTotalArray{"$SecondLook"}};
	}
#end of little sub to get $SecondLook
#########################################################


#########################################################
# Sub to look for find assumptions data
sub AssumptSub {
	#initialise
	open(ASSUMPT, '<', $Assumptfile);
	$FileToClean = "ASSUMPT";
	&ImportFile;
	@AssumptNEXUSstore = @OriginalNEXUSstore;
	@Assumptions = @FormatNEXUSstore;
	$loopend = 0;
	@Typesets =();
	@CodonSets = ();
	$AreThereCodons = 0;

	#Assess characters
	$TypesetCheck = 0;
	$CTypeCheck = 0;

	for($i = 0; $i <= $#AssumptNEXUSstore; $i++){
		$TempIn = $AssumptNEXUSstore[$i];
		chomp($TempIn);
		if($TempIn =~ /^typeset$/i){
			$i++;
			$TempIn = $AssumptNEXUSstore[$i];
			chomp($TempIn);

			if($TempIn eq "*"){
				$i++;
				$TempIn = $TempIn . $AssumptNEXUSstore[$i];
				chomp($TempIn);
				}

			push(@Typesets, $TempIn);
			}
		if($TempIn =~ /^ctype$/i){
			$CTypeCheck++;
			}
		if($DataType =~ /dna|rna|nucleotide/i){
			if($TempIn =~ /^codonposset$/i){
				$i++;
				$TempIn = $AssumptNEXUSstore[$i];
				chomp($TempIn);

				if($TempIn eq "*"){
					$i++;
					$TempIn = $AssumptNEXUSstore[$i];
					chomp($TempIn);
					}
				push(@CodonSets, $TempIn);
				$AreThereCodons = 1;
				}
			}
		}
	}
#end of AssumptSub
#########################################################


#########################################################
#sub to check to see if there are any default character types
sub DefAssumptSub {
	@DeftypeList = ();
	for($i = 0; $i < $#Assumptions; $i++){
		$typetext = $Assumptions[$i];
		if($typetext =~ /deftype/i){
			$DefTypeFromNexus = $Assumptions[$i + 2];
			$DefTypeFromNexus =~ tr/A-Z/a-z/;
			if($DefTypeFromNexus =~ /dollo$/){$DefTypeFromNexus = "dollo.up"}
			if($DefTypeFromNexus =~ /irrev$/){$DefTypeFromNexus = "irrev.up"}
			push(@DeftypeList, $DefTypeFromNexus);
			}
		}

	if(@DeftypeList > 1){
		$SetType = $DeftypeList[$#DeftypeList];
		}
	elsif(@DeftypeList == 1){
		$SetType = $DeftypeList[$#DeftypeList];
		}
	else{
		$SetType = "unord";
		}
	}
#end of sub to check to see if there are any default character types
#########################################################


#########################################################
#Sub to set default types
sub SettingType {
	@CharTypes2 = ();
	@Understood = ();
	for($i = 0; $i < $NCHAR; $i++){$CharTypes2[$i] = $SetType};
	for($i = 0; $i < $NCHAR; $i++){$Understood[$i] = $Misunderstood}
	}
#end of sub to set default types
#########################################################


#########################################################
#sub to examine and read in typesets
sub ExaminationTypeset {
	#Big While loop to examine Typesets
	while(1){# Start of Big While loop to get Typesets
		if(@Typesets == 0){
			$Insert = "No Typesets detected - n";
			$Insert2 =  "No Typesets detected - ";
			&NoTypeSets ;
			last
			}

		#BIG LOOP if there are typesets...

		#Pick out default typesets (marked with a *)
		else{&TypeSetSearch}

		#what to do if there are more than one typeset...
		if($TypesetCheck == 1){
			while(1){
				$Default  = @DefaultTypesets;
				$Other = @OtherTypesets;
				if(($Default == 0) and ($Other > 1)){
					print "\nNo default Typeset detected\nPlease select a typeset from this list:\n";
					$TypeNumber = 0;
					@TestList = ();
					foreach $i(@OtherTypesets){
						$TypeNumber ++;
						push(@TestList, " $TypeNumber ");
						print "   $TypeNumber. $i\n"
						};
					$TypeNumber ++;
					print "   $TypeNumber. Set a global type \(e.g. All ordered, All unordered, etc\)\n";
					print "      or use a Nexus format assumption file\n";
					push(@TestList, " $TypeNumber ");
					print "\n   ";
					$TypeNumber2 = <stdin>;
					chomp ($TypeNumber2);
					$TestInput = grep(/ $TypeNumber2 /, @TestList);
					if(($TestInput == 0) or (length($TypeNumber2) == 0)){next}
					elsif($TypeNumber2 == $TypeNumber){
						$DefNoDef = 2;
						last
						}
					else {$DefNoDef = 0; last}
					}

				elsif($Default == 1){
					print "\nOne default Typeset detected \($DefaultTypesets[0]\)\nDo you wish to use it?:\n";
					$TypeNumber = 0;
					@TestList = ();
					$TypeNumber ++;
					push(@TestList, " $TypeNumber ");
					print "   $TypeNumber. $DefaultTypesets[0]\n";
					$TypeNumber ++;
					print "   $TypeNumber. Set a global type (e.g. All ordered, All unordered, etc\)\n";
					print "      or use a Nexus format assumption file\n";
					push(@TestList, " $TypeNumber ");
					print "\n   ";
					$TypeNumber2 = <stdin>;
					chomp ($TypeNumber2);
					$TestInput = grep(/ $TypeNumber2 /, @TestList);
					if(($TestInput == 0) or (length($TypeNumber2) == 0)){next}
					elsif($TypeNumber2 == $TypeNumber){
						$DefNoDef = 2;
						last
						}
					else {$DefNoDef = 1; last}
					}

				elsif($Default > 1){
					print "\nMore than one default Typeset detected (i.e. name preceeded by an *)\nPlease select one from this list:\n";
					$TypeNumber = 0;
					@TestList = ();
					foreach $i(@DefaultTypesets){
						$TypeNumber ++;
						push(@TestList, $TypeNumber);
						print "   $TypeNumber. $i\n";
						};
					$TypeNumber ++;
					print "   $TypeNumber. Set a global type (e.g. All ordered, All unordered, etc\)\n";
					print "      or use a Nexus format assumption file\n";					push(@TestList, $TypeNumber);
					print "\n   ";
					$TypeNumber2 = <stdin>;
					chomp ($TypeNumber2);
					$TestInput = grep(/$TypeNumber2/, @TestList);

					if(($TestInput == 0) or (length($TypeNumber2) == 0)){next}
					elsif($TypeNumber2 == $TypeNumber){
						$DefNoDef = 2;
						last
						}
					else {$DefNoDef = 1; last}
					}

				else{
					print "\nOne typeset detected \($OtherTypesets[0]\)\nDo you wish to use it?:\n";
					$TypeNumber = 0;
					@TestList = ();
					$TypeNumber ++;
					push(@TestList, " $TypeNumber ");
					print "   $TypeNumber. $OtherTypesets[0]\n";
					$TypeNumber ++;
					print "   $TypeNumber. Set a global type (e.g. All ordered, All unordered, etc\)\n";
					print "      or use a Nexus format assumption file\n";
					push(@TestList, " $TypeNumber ");
					print "\n   ";
					$TypeNumber2 = <stdin>;
					chomp ($TypeNumber2);
					$TestInput = grep(/ $TypeNumber2 /, @TestList);
					if(($TestInput == 0) or (length($TypeNumber2) == 0)){next}
					elsif($TypeNumber2 == $TypeNumber){
						$DefNoDef = 2;
						last
						}
					else {$DefNoDef = 0; last}
					}

				}#end of if loop
			}#end of while loop

	#Pick Typeset to use
		if($DefNoDef == 0){
			$TypeNumber2 --;
			$UsedTypeset = $OtherTypesets[$TypeNumber2]
			}
		elsif($DefNoDef == 1){
			$TypeNumber2 --;
			$UsedTypeset = $DefaultTypesets[$TypeNumber2]
			}
		elsif($DefNoDef == 2){
			$Insert = "User rejected typesets\n\nN";
			$Insert2 =  "";
			&NoTypeSets ;
			last
			}

		$Quotey = "";
		unless($UsedTypeset =~ /^\'|\'$/){
			$Quotey = "\'";
			}

		print OUTPUT "\n$Quotey<b>$UsedTypeset</b>$Quotey used as typeset:<br><br>\n";

	#####################################################
	# Move to selected typeset
		#find which number type set it is in the NEXUS
		$WhichTypeSet = 0;
		for($t = 0; $t <= $#Typesets; $t++){
			$WhichTypeSet++;
			if($Typesets[$t] eq $UsedTypeset){
				last;
				}
			}

		$i = 0;
		$TypeCounter = 0;
		while($i <= $#AssumptNEXUSstore){
			$TypeSetInput = $AssumptNEXUSstore[$i];
			chomp($TypeSetInput);
			if($TypeSetInput =~ /^typeset$/i){
				$TypeCounter++;
				if($TypeCounter == $WhichTypeSet){
					until($TypeSetInput eq '='){
						$i++;
						$TypeSetInput = $AssumptNEXUSstore[$i];
						chomp($TypeSetInput);
						}
					$i--;
					$LineIn = $AssumptNEXUSstore[$i];
					last;
					}
				}
			$i++;
			}

	#####################################################
	# Start to read in character types
		$TypeTest = 0; #lets it look for type

		#if vector input
		if($AssumptNEXUSstore[$i] =~ /^\(vector\)$/i){
			$i += 2;
			$LineIn = $AssumptNEXUSstore[$i];

			$CharCount = 0;
			until($LineIn =~ /;/){
				$LineIn =~ tr/A-Z/a-z/;
				chomp($LineIn);
				$FinalType = &SetTypeSub;

				$Understood[$CharCount] =  $Misunderstood;
				$CharTypes2[$CharCount] = $FinalType;
				$i++;
				$CharCount++;
				$LineIn = $AssumptNEXUSstore[$i];
				}
			}

		#if standard input
		elsif($AssumptNEXUSstore[$i] =~ /^\(standard\)$/i){
			$i += 2;
			&ReadInTypes2;
			}
		else{
			$i += 2;
			&ReadInTypes2;
			}
		last
		}#End of Big While loop to get Typesets
	}
#end of sub ExaminationTypeset
#########################################################


#########################################################
#Sub for reading codon blocks
sub CodonReader {
	$CodonTypes = 0;

	#Big While loop to examine Codon sets
	while(1){# Start of Big While loop to get Codon sets
		if(@CodonSets == 0){
			last;
			}

		#BIG LOOP if there are Codon sets...

		#what to do if there are more than one typeset...
		if(@CodonSets >= 1){
			while(1){
				$Default = @CodonSets;

				if($Default == 1){
					$CodonTypes = 1;
					last;
					}

				elsif($Default > 1){
					print "\nMore than one CodonPosSet detected\nPlease select one from this list:\n";
					$TypeNumber = 0;
					@TestList = ();
					foreach $i(@CodonSets){
						$TypeNumber ++;
						push(@TestList, $TypeNumber);
						print "   $TypeNumber. $i\n";
						};
					$TypeNumber2 = <stdin>;
					chomp ($TypeNumber2);
					$TestInput = grep(/$TypeNumber2/, @TestList);

					if(($TestInput == 0) or (length($TypeNumber2) == 0)){next}
					else {
						$CodonTypes = 1;
						last;
						}
					}
				}#end of while loop
			}#end of if loop


		#Pick Typeset to use
		$TypeNumber2 --;
		$UsedCodonSet = $CodonSets[$TypeNumber2];

		#make default setting
		@UsedCodons = ();
		for($Z = 0; $Z < $NCHAR; $Z++){
			push(@UsedCodons, '?');
			}

		#####################################################
		# Move to selected codon set
		#find which number type set it is in the NEXUS
		$WhichCodonSet = 0;
		for($t = 0; $t <= $#CodonSets; $t++){
			$WhichCodonSet++;
			if($CodonSets[$t] eq $UsedCodonSet){
				last;
				}
			}
		$i = 0;
		$TypeCounter = 0;
		while($i <= $#AssumptNEXUSstore){
			$CodonSetInput = $AssumptNEXUSstore[$i];
			chomp($CodonSetInput);
			if($CodonSetInput =~ /^codonposset$/i){
				$TypeCounter++;
				if($TypeCounter == $WhichCodonSet){
					until($CodonSetInput eq '='){
						$i++;
						$CodonSetInput = $AssumptNEXUSstore[$i];
						chomp($CodonSetInput);
						}
					$i--;
					$LineIn = $AssumptNEXUSstore[$i];
					last;
					}
				}
			$i++;
			}


		#####################################################
		# Start to read in character types
		$TypeTest = 0; #lets it look for type

		#if vector input
		if($AssumptNEXUSstore[$i - 1] =~ /^\(vector\)$/i){
			$i++;
			$LineIn = $AssumptNEXUSstore[$i];

			$CharCount = 0;
			until($LineIn =~ /;/){
				$LineIn =~ tr/A-Z/a-z/;
				chomp($LineIn);
				$FinalType = &SetTypeSub;
				$UsedCodons[$CharCount] = $FinalType;
				$i++;
				$CharCount++;
				$LineIn = $AssumptNEXUSstore[$i];
				}
			}

		#if standard input
		elsif($AssumptNEXUSstore[$i] =~ /^\(standard\)$/i){
			$i += 2;
			&ReadInTypes2;
			}
		else{
			$i += 2;;
			&ReadInTypes2;
			}

		@Tempy = ();
		for($Ru = 0; $Ru < $NCHAR; $Ru++){
			push(@Tempy, '0');
			}
		@tempN = ('charset', @Tempy);
		@temp1 = ('charset', @Tempy);
		@temp2 = ('charset', @Tempy);
		@temp3 = ('charset', @Tempy);
		@tempWhat = ('charset', @Tempy);

		for($FG = 1; $FG <= $NCHAR; $FG++){
			if($UsedCodons[$FG - 1] eq 'n'){
				$tempN[$FG] = 1;
				}
			if($UsedCodons[$FG - 1] eq '1'){
				$temp1[$FG] = 1;
				}
			if($UsedCodons[$FG - 1] eq '2'){
				$temp2[$FG] = 1;
				}
			if($UsedCodons[$FG - 1] eq '3'){
				$temp3[$FG] = 1;
				}
			if($UsedCodons[$FG - 1] eq '?'){
				$tempWhat[$FG] = 1;
				}
			}

		$TotalCharLabels{"noncoding"} = [@tempN];
		$TotalCharLabels{"pos1"} = [@temp1];
		$TotalCharLabels{"pos2"} = [@temp2];
		$TotalCharLabels{"pos3"} = [@temp3];
		$TotalCharLabels{"?"} = [@tempWhat];

		$CodonTypes = 0;
		last
		}#end of big while loop to get codons
	}
#end of sub for reading codon blocks
#########################################################


#########################################################
#Subroutine for getting the type of character
sub SetTypeSub {

	$LineIn =~ tr/A-Z/a-z/;
	if($LineIn =~ /:$/){chop($LineIn)}

	if($LineIn eq "ord"){
		$Misunderstood = 0;
		return("ord")
		}
	elsif($LineIn eq "unord"){
		$Misunderstood = 0;
		return("unord")
		}
	elsif(($LineIn eq "irrev") or ($LineIn eq "irrev.up")){
		$Misunderstood = 0;
		return("irrev.up")
		}
	elsif($LineIn eq "irrev.dn"){
		$Misunderstood = 0;
		return("irrev.dn")
		}
	elsif(($LineIn eq "dollo") or ($LineIn eq "dollo.up")){
		$Misunderstood = 0;
		return("dollo.up")
		}
	elsif($LineIn eq "dollo.dn"){
		$Misunderstood = 0;
		return("dollo.dn")
		}
	elsif($CodonTypes == 1){
		if($LineIn eq "n"){
			$Misunderstood = 0;
			return("n")
			}
		elsif($LineIn eq "1"){
			$Misunderstood = 0;
			return("1")
			}
		elsif($LineIn eq "2"){
			$Misunderstood = 0;
			return("2")
			}
		elsif($LineIn eq "3"){
			$Misunderstood = 0;
			return("3")
			}
		elsif($LineIn eq "?"){
			$Misunderstood = 0;
			return("?")
			}
		}
	else{
		$Misunderstood = 1;
		return("unord")
		}
	};
#end of sub for getting the type of character
#########################################################


#########################################################
#sub to sort names of typesets in default and other
sub TypeSetSearch {
 	$TypesetCheck = 1;
 	@DefaultTypesets = ();
 	@OtherTypesets = ();
 	foreach $i(@Typesets){
 		if($i =~ /^\*/){push (@DefaultTypesets, $i)}
 		else{push (@OtherTypesets, $i)}
 		}
	}
#end of sub to sort names of typesets in default and other
#########################################################


#########################################################
#sub to read in typesets
sub ReadInTypes2 {
	$IncreNot = 0;

	#set remainder
	@TipTip = ('charset');
	for($yt = 0; $yt <= $NCHAR; $yt++){
		push(@TipTip, 1);
		}
	$TotalCharLabels{"remainder"} = [@TipTip];

	#big loop to read types
	until($LineIn =~ /;/){
		$LineIn = $AssumptNEXUSstore[$i];
		chomp($LineIn);
		$LineIn =~ tr/A-Z/a-z/;

		if($LineIn =~ /^remainder$/i){$LineIn = "remainder"}
		if($LineIn =~ /^all$/i){$LineIn = "all"}
		if($LineIn =~ /^uninf$/i){$LineIn = "uninf"}
		if($LineIn =~ /^constant$/i){$LineIn = "constant"}
		if($LineIn =~ /^gapped$/i){$LineIn = "gapped"}
		if($LineIn =~ /^missambig$/i){$LineIn = "missambig"}
		@FirstOne = @{$TotalCharLabels{"$LineIn"}};

		$NextLine = $AssumptNEXUSstore[$i + 1];
		chomp($NextLine);
		$NextLine =~ tr/A-Z/a-z/;
		if($NextLine =~ /^remainder$/i){$NextLine = "remainder"}
		if($NextLine =~ /^all$/i){$NextLine = "all"}
		if($NextLine =~ /^uninf$/i){$NextLine = "uninf"}
		if($NextLine =~ /^constant$/i){$NextLine = "constant"}
		if($NextLine =~ /^gapped$/i){$NextLine = "gapped"}
		if($NextLine =~ /^missambig$/i){$NextLine = "missambig"}

		#Look for type
		if($TypeTest == 0){
			$FinalType = &SetTypeSub; #Go to sub which reads type
			$TypeTest = 1;
			$i++;
			$QuickTest = $AssumptNEXUSstore[$i];
			chomp($QuickTest);
			if($QuickTest eq ":"){
				$i++;
				}
			next
			}

		#Look for characters  $Misunderstood2 = 0; @CharTypes2
		if($TypeTest == 1){
			#If "all"
			if($LineIn =~ /^all$/i){
				for($y = 1; $y <= $NCHAR; $y++){
					if($CodonTypes == 1){$UsedCodons[$y - 1] = $FinalType;}
					else{
						$CharTypes2[$y - 1] = $FinalType;
						$Understood[$y - 1] = $Misunderstood;
						}
					$TotalCharLabels{"remainder"}[$y] = 0;
					$CTypeChanges = 1;
					}
				last
				}

			#after incremented span
			elsif($LineIn =~ /,/){
				$i++;
				$TypeTest = 0;
				next
				}

			#if number at end of one type
			elsif(($FirstOne[0] eq "char") and ($NextLine =~ /,/)){
				if($CodonTypes == 1){$UsedCodons[$FirstOne[1] - 1] =  $FinalType}
				else{
					$CharTypes2[$FirstOne[1] - 1] =  $FinalType;
					$CTypeChanges = 1;
					$Understood[$FirstOne[1] - 1] = $Misunderstood;
					}
				$TotalCharLabels{"remainder"}[$FirstOne[1]] = 0;
				$i = $i + 2;
				$TypeTest = 0;
				next
				}

			#if charset end of one type
			elsif(($FirstOne[0] eq "charset") and ($NextLine =~ /,/)){
				shift(@FirstOne);
				$County = 1;
				foreach $qr (@FirstOne){
					if($qr == 1){
						if($CodonTypes == 1){$UsedCodons[$County - 1] =  $FinalType}
						else{
							$CharTypes2[$County - 1] =  $FinalType;
							$CTypeChanges = 1;
							$Understood[$County - 1] = $Misunderstood;
							}
						$TotalCharLabels{"remainder"}[$County] = 0;
						}
					$County++;
					}
				$i = $i + 2;
				$TypeTest = 0;
				next
				}

			#if a span of characters
			elsif($NextLine eq "-"){
				$NextNextLine = $AssumptNEXUSstore[$i + 2];
				@ThirdOne = @{$TotalCharLabels{"$NextNextLine"}};

				$IncreNot = 0;
				$Incrementer = $AssumptNEXUSstore[$i + 3];
				if($Incrementer =~ /^\\$/){
					$Incrementer = $AssumptNEXUSstore[$i + 4];
					if($Incrementer =~ /^\d$/){
						$IncreNot = 1;
						}
					else{
						$Incrementer = 1;
						$IncreNot = 0;
						}
					}
				else{$Incrementer = 1;}


				#if all
				if($NextNextLine eq "."){
					for($z = $FirstOne[1]; $z <= $NCHAR; $z += $Incrementer){
						if($CodonTypes == 1){$UsedCodons[$z - 1] =  $FinalType}
						else{
							$CharTypes2[$z - 1] =  $FinalType;
							$CTypeChanges = 1;
							$Understood[$z - 1] = $Misunderstood;
							}
						$TotalCharLabels{"remainder"}[$z] = 0;
						}
					}

				else{
					for($z = $FirstOne[1]; $z <= $ThirdOne[1]; $z += $Incrementer){
						if($CodonTypes == 1){$UsedCodons[$z - 1] =  $FinalType}
						else{
							$CharTypes2[$z - 1] =  $FinalType;
							$CTypeChanges = 1;
							$Understood[$z - 1] = $Misunderstood;
							}
						$TotalCharLabels{"remainder"}[$z] = 0;
						}
					}
				if($IncreNot == 1){$i += 4}
				$i++;
				next
				}

			#if just a number
			elsif(($FirstOne[0] eq "char") and ($NextLine ne "-")){
				#Remainder
				if($CodonTypes == 1){$UsedCodons[$FirstOne[1] - 1] =  $FinalType}
				else{
					$CharTypes2[$FirstOne[1] - 1] =  $FinalType;
					$CTypeChanges = 1;
					$Understood[$FirstOne[1] - 1] =  $Misunderstood;
					}
				$TotalCharLabels{"remainder"}[$FirstOne[1]] = 0;
				$i++;
				next
				}

			elsif(($FirstOne[0] eq "charset") and ($NextLine ne "-")){
				shift(@FirstOne);
				$County = 1;
				foreach $qr (@FirstOne){
					if($qr == 1){
						if($CodonTypes == 1){$UsedCodons[$County - 1] =  $FinalType}
						else{
							$CharTypes2[$County - 1] =  $FinalType;
							$CTypeChanges = 1;
							$Understood[$County - 1] = $Misunderstood;
							}
						$TotalCharLabels{"remainder"}[$County] = 0;
						}
					$County++;
					}
				$i++;
				next
				}

			#in case of a space or whatever
			else{
				$i++;
				next
				}
			}#end of looking for characters
		}#end of big typeset loop!
	}
#end of sub ReadInType!
#########################################################


#########################################################
#If there are no type sets...
sub NoTypeSets {
	while(1){#$SetType"
		print "\n$Insert";
		print "ame a Nexus format assumptions file, or assume\n";

		print "   1 All unordered\n";
		print "   2 All ordered\n";
		print "   3 All Dollo up\n";
		print "   4 All Dollo down\n";
		print "   5 All irreverable up\n";
		print "   6 All irreverable down\n";
		print "   (current default is \'$SetType\')\n   ";

		$Answer = <stdin>;
		chomp ($Answer);

		if($Answer =~ /^q$/){die}
		elsif($Answer == 1){$SetType = "unord";
			$fullset = "unordered";
			&SettingType;
			$LoopBreak = 1
			}
		elsif($Answer == 2){$SetType = "ord";
			$fullset = "ordered";
			&SettingType;
			$LoopBreak = 1
			}
		elsif($Answer == 5){
			$SetType = "irrev.up";
			$fullset = "irreverable (up)";
			&SettingType;
			$LoopBreak = 1
			}
		elsif($Answer == 6){
			$SetType = "irrev.dn";
			$fullset = "irreverable (down)";
			&SettingType;
			$LoopBreak = 1
			}
		elsif($Answer == 3){
			$SetType = "dollo.up";
			$fullset = "Dollo (up)";
			&SettingType;
			$LoopBreak = 1
			}
		elsif($Answer == 4){
			$SetType = "dollo.dn";
			$fullset = "Dollo (down)";
			&SettingType;
			$LoopBreak = 1
			}
		else{
			$Assumptfile = $Answer;
			if(-e $Assumptfile){
				&AssumptSub;
				$SetType = "unord";
				$Misunderstood = 0;
				&SettingType;
				&ExaminationTypeset;
				$LoopBreak = 1;
				last
				}
			else{
				print "\nThe file \'$Assumptfile\' could not be found.  Please try again\n";
				$Insert = "N";
				next
				}
			}

		print OUTPUT "\n<p>" . $Insert2 . "User set all characters to be $fullset\<br><br></p>\n";
		$UsedTypeset = "User_Set";
		last
		}#end of input loop

	if($LoopBreak == 1){last}
	else{next}
	}

#end of sub NoTypeSets
#########################################################


#########################################################
#sub to input ctype statements
sub CtypeReadin {
	$CTypeChanges = 0;
	for($i = 0; $i <= $#AssumptNEXUSstore; $i++){
		$LineIn = $AssumptNEXUSstore[$i];
		chomp($LineIn);
		#if you find a ctype, read it in
		if($LineIn =~ /^ctype$/i){
			$i++;
			$TypeTest = 0;
			&ReadInTypes2;
			}
		}
	#if ctype has changed @chartype2, it overrides any charsets which may be present
	if($CTypeChanges == 1){
		$UsedTypeset = "\'User Set\'";
		$TwoForTea = "";
		$DoubleTrouble = "a ";
		if($CTypeCheck > 1){
			$TwoForTea = "s";
			$DoubleTrouble = "";
			}

		print OUTPUT "\n<P>User set types with $DoubleTrouble\'CType\' command$TwoForTea:<br><br></p>\n";
		}
	}
#end ofsub to input ctype statements
#########################################################


#########################################################
#sub to find exsets and exclude characters
sub ExSetFinder {
	@ExSetList = ();
	$JJ = 0;
	for($EE = 0; $EE <= $#ExsetNexus; $EE++){
		$yy = $ExsetNexus[$EE];
		$JJ = $ExsetNexus[$EE + 1];
		chomp($yy);
		chomp($JJ);

		if($yy =~ /^exset$/i){
			if($JJ eq "*"){
				$EE++;
				$JJ = $ExsetNexus[$EE + 1];
				}
			push(@ExSetList, $JJ);
			}
		}

	#Remove asterisk from set names...
	for($RT = 0; $RT <= $#ExSetList; $RT++){
		$ExSetList[$RT] =~ s/\*//g;
		}

	#####################################################
	# Start of Big While loop to get Exsets
	$AllIncTest = 0;
	while(1){
		if(@ExSetList == 0){last}
		elsif(@ExSetList > 1){&ExSetSearch}
		else{
			while(1){
				print "\nOne ExSet detected \($ExSetList[0]\).  Do you want to use it?";
				print "\n1. $ExSetList[0]";
				print "\n2. All characters to be included \(overrides \'eliminate\' statements\)\n   ";

				$All_in_or_not = <STDIN>;
				chomp $All_in_or_not;

				if($All_in_or_not == 1){last}
				elsif($All_in_or_not == 2){last}

				}
			if($All_in_or_not == 1){
				$WhichExsetNumber = 1;
				$UsedExset = $ExSetList[0];
				$WhereToGo = 1;
				}
			if($All_in_or_not == 2){$AllIncTest = 1; last}
			}#goto read-in

		if($WhereToGo == 1){
			&GotoExSet
			};

		#Remove any spurious Exset codings
		for($DD = 0; $DD <= $#FinalExSet; $DD++){
			if(($FinalExSet[$DD] != 0) and ($FinalExSet[$DD] != 1)){
				print "\rExset data for Character ";
				print $DD + 1;
				print " is not understood.  Character excluded\r";
				$FinalExSet[$DD] = 1
				}
			}
		last
		}#End of Big While loop to get Exsets
		#####################################################
		#####################################################

	#look for eliminate statements
	&EliminateStatement;

	##############################
	#read list and Exclude characters...


	for($GT = 0; $GT <= $#NexusArray; $GT++){
		for($ED = 1; $ED <= $NCHAR; $ED++){
			if($FinalExSet[$ED-1] == 1){
				$NexusArray[$GT][$ED] = $Missing
				}
			}
		}
	}
#end of sub to find exsets and exclude characters
#########################################################


#########################################################
#sub to find & implememt 'eliminate' statement...
sub EliminateStatement {
	unless($AllIncTest == 1){
		@EliminateList = ();
		@EliminateDiff = ();
		for($j = 0; $j < $NCHAR; $j++){
			$EliminateList[$j] = 0;
			}

		for($i = 0; $i <= $#ExsetNexus; $i++){
			$LineIn = $ExsetNexus[$i];
			if($LineIn =~ /^eliminate$/i){
				$i++;
				until($LineIn eq ";"){

					&ExsetNextLine;

					$ik = $i + 1;
					&ExsetNextLineTwo;

					#if a charset...
					if($WhatChar[0] eq "charset"){
						shift(@WhatChar);
						$wf = 0;
						foreach $wd (@WhatChar){
							if($wd == 1){
								$EliminateList[$wf] = 1;
								}
							$wf++;
							}
						$i++;
						next;
						}

					elsif($WhatChar[0] eq "char"){
						#number and number
						if($NextLine ne "-"){
							$EliminateList[$WhatChar[1] - 1] = 1;
							$i++;
							next;
							}

						#number and dash
						elsif($NextLine eq "-"){
							$ik++;
							&ExsetNextLineTwo;
							$i = $ik+1;

							$IncreNot = 0;
							$Incrementer = $ExsetNexus[$i];
							if($Incrementer =~ /^\\$/){
								$Incrementer = $ExsetNexus[$i + 1];
								if($Incrementer =~ /^\d$/){
									$IncreNot = 1;
									}
								else{
									$Incrementer = 1;
									$IncreNot = 0;
									}
								}
							else{$Incrementer = 1;}

							#dash and number
							if($WhatCharTwo[0] eq "char"){
								for($k = $WhatChar[1]; $k <= $WhatCharTwo[1]; $k += $Incrementer){
									$EliminateList[$k -1] = 1;
									}
								}

							if($NextLine eq "."){
								for($k = $WhatChar[1]; $k <= $NCHAR; $k += $Incrementer){
									$EliminateList[$k -1] = 1;
									}
								}
							if($IncreNot == 1){$i += 2}
							next;
							}
						}
					else{next}
					}
				}
			}
		for($L = 0; $L <= $#EliminateList; $L++){
			$InOut = $EliminateList[$L];
			if($InOut == 1){
				$SomeEliminated = 1;
				$ExsetInOut = $FinalExSet[$L];
				if($ExsetInOut == 0){
					$FinalExSet[$L] = 1;
					push(@EliminateDiff, ($L + 1));
					}
				}
			}
		}
 	}

#end of sub to find & implememt 'eliminate' statement...
#########################################################


#########################################################
#Reads in standard-format Exsets
sub ReadInTypes {
	$IncreNot = 0;
	#big loop to read types
	until($LineIn eq "\;"){

		&ExsetNextLine;
		$ik = $i+1;
		&ExsetNextLineTwo;

		#if a charset
		if($WhatChar[0] eq "charset"){
			shift(@WhatChar);
			$tg = 0;
			foreach $fr (@WhatChar){
				if($fr == 1){
					$FinalExSet[$tg] = 1;
					}
				$tg++;
				}
			$i++;
			next;
			}

		elsif($WhatChar[0] eq "char"){
			#if a span of characters
			if($NextLine eq "-"){
				$ik++;
				&ExsetNextLineTwo;

				$IncreNot = 0;
				$Incrementer = $ExsetNexus[$ik + 1];
				if($Incrementer =~ /^\\$/){
					$Incrementer = $ExsetNexus[$ik + 2];
					if($Incrementer =~ /^\d$/){
						$IncreNot = 1;
						}
					else{
						$Incrementer = 1;
						$IncreNot = 0;
						}
					}
				else{$Incrementer = 1;}

				# if span is till end
				if($NextLine eq "."){
					for($z = $WhatChar[1]; $z <= $NCHAR; $z += $Incrementer){
						$FinalExSet[$z - 1] =  1;
						}
					}

				# if limited span
				else{
					for($z = $WhatChar[1]; $z <= $WhatCharTwo[1]; $z += $Incrementer){
						$FinalExSet[$z - 1] =  1;
						}
					}
				$i = $ik + 1;
				if($IncreNot == 1){$i += 2}
				next
				}

			#if just a char
			else{
				$FinalExSet[$WhatChar[1] - 1] =  1;
				$i++;
				next
				}
			}

		#terminate
		elsif($LineIn eq "\;"){last}

		}#end of looking for characters
	}
#end sub to read standard-format exsets
#########################################################


#########################################################
#Sub if there several are Exsets...
sub ExSetSearch {
	while(1){
		print "\nMore than one Exset detected\nPlease select an Exset from this list:\n";
		$WhereToGo =0;
		$TypeNumber = 0;
		@TempList = ();
		foreach $i(@ExSetList){
			$TypeNumber ++;
			push(@TempList, $TypeNumber);
			print "   $TypeNumber. $i\n"
			};
		$TypeNumber ++;
		print "   $TypeNumber. All characters to be included \(overrides \'eliminate\' statements\)\n";
		push(@TempList, $TypeNumber);
		print "\n   ";
		$TypeNumber3 = <stdin>;
		chomp($TypeNumber3);
		$TestInput = grep(/$TypeNumber3/, @TempList);

		if(($TestInput == 0) or (length($TypeNumber3) == 0)){next}
		#if all inc.
		elsif($TypeNumber3 == $TypeNumber){
			$AllIncTest = 1;
			last;
			}
		#if exsetnumber
		else{
			$WhichExsetNumber = $TypeNumber3;
			$UsedExset = $ExSetList[$TypeNumber3 - 1];
			$WhereToGo = 1;
			last
			}
		}
	}#end of sub if there several are Exsets...
#########################################################


#########################################################
# Locate Exset & read in data
sub GotoExSet {
	# Move to selected typeset
	$FindCheck = 0;
	for($i = 0; $i <= $#ExsetNexus; $i++){
		$TestLine = $ExsetNexus[$i];
		chomp($TestLine);

		if($TestLine =~ /^exset$/i){
			$FindCheck++;
			}
		if($FindCheck == $WhichExsetNumber){
			until($TestLine eq "="){
				$i++;
				$TestLine = $ExsetNexus[$i];
				chomp($TestLine);
				}
			$i--;
			$TestLine = $ExsetNexus[$i];
			chomp($TestLine);
			$FindCheck = 2;
			last;
			}
		}

	#####################################################
	# Start to read in character types

	#Reads in vector-format Exsets
	if($ExsetNexus[$i] =~ /^\(vector\)$/i){

		@FinalExSet = ();
		$CharCount = 0;
		$i += 2;
		until($CharCount == $NCHAR){
			$LineIn = $ExsetNexus[$i];
			chomp($LineIn);
			@TempExSet = ();
			until(length($LineIn) == 0){
				$ExIn = chop($LineIn);
				unshift(@TempExSet, $ExIn);
				$CharCount++;
				}
			push(@FinalExSet,@TempExSet);
			$i ++;
			}
		}

	#Identifies standard-format Exsets
	elsif($ExsetNexus[$i] =~ /^\(standard\)$/i){
		$i += 2;
		$LineIn = $ExsetNexus[$i];
		&ReadInTypes;
		}
	else{
		$i +=2;
		$LineIn = $ExsetNexus[$i];
		&ReadInTypes;
		}
	}
#end of sub GotoExSet
#########################################################


#########################################################
#little sub to get next Exset line
sub ExsetNextLine {
	$LineIn = $ExsetNexus[$i];
	chomp($LineIn);

	if($LineIn =~ /^all$/i){
		$LineIn = "all";
		}
	if($LineIn =~ /^uninf$/i){
		$LineIn = "uninf";
		}
	if($LineIn =~ /^constant$/i){
		$LineIn = "constant";
		}
	if($LineIn =~ /^gapped$/i){
		$LineIn = "gapped";
		}
	if($LineIn =~ /^missambig$/i){
		$LineIn = "missambig";
		}
	if($LineIn =~ /^pos1$/i){
		$LineIn = "pos1";
		}
	if($LineIn =~ /^pos2$/i){
		$LineIn = "pos2";
		}

	if($LineIn =~ /^pos3$/i){
		$LineIn = "pos3";
		}
	if($LineIn =~ /^noncoding$/i){
		$LineIn = "noncoding";
		}
	@WhatChar = @{$TotalCharLabels{"$LineIn"}};

	}
#little sub to get next Exset line
#########################################################


#########################################################
#little sub to get next Exset line but one
sub ExsetNextLineTwo {
	$NextLine = $ExsetNexus[$ik];
	chomp($NextLine);

	if($NextLine =~ /^all$/i){
		$NextLine = "all";
		}
	if($NextLine =~ /^uninf$/i){
		$NextLine = "uninf";
		}
	if($NextLine =~ /^constant$/i){
		$NextLine = "constant";
		}
	if($NextLine =~ /^gapped$/i){
		$NextLine = "gapped";
		}
	if($NextLine =~ /^missambig$/i){
		$NextLine = "missambig";
		}
	if($LineIn =~ /^pos1$/i){
		$LineIn = "pos1";
		}
	if($LineIn =~ /^pos2$/i){
		$LineIn = "pos2";
		}
	if($LineIn =~ /^pos3$/i){
		$LineIn = "pos3";
		}
	if($LineIn =~ /^noncoding$/i){
		$LineIn = "noncoding";
		}
	@WhatCharTwo = @{$TotalCharLabels{"$NextLine"}};
	}
#little sub to get next Exset line but one
#########################################################


#########################################################
#sub to get the NEXUS file input in basic form
sub ImportFile {
	#Put whole infile into a list
	@InputNexus = ();
	@MacPCNexus = ();
	@InputNexus = <$FileToClean>;

	#SET LINE BREAKS TO LOCAL TYPE
	$Mac = chr(015);
	$PC = chr(012);
	for ($i = 0; $i <= $#InputNexus; $i++) {
		$line = $InputNexus[$i];
		$line =~ s/$Mac|$PC/\n/g;
		$line =~ tr/\n/\n/s;
		@TempList = split(/\n/, $line);
		unless(@TempList == 0){
			for($NT = 0; $NT <= $#TempList; $NT++){
				$TempList[$NT] .= "\n";
				}

			push(@MacPCNexus, @TempList);
			}
		}

	@FormatNEXUSstore = ();
	@NewBracketNEXUSstore = ();
	@OriginalNEXUSstore = ();
	@BracketNEXUSstore = ();

	#Set spaces to single
	for ($i = 0; $i <= $#MacPCNexus; $i++) {
		$TheLine = $MacPCNexus[$i];

		$TheLine =~ s/\\/ \\ /g;
		$TheLine =~ s/\// \/ /g;
		$TheLine =~ s/\[/ [ /g;
		$TheLine =~ s/\]/ ] /g;
		$TheLine =~ s/;/ ; /g;
		$TheLine =~ s/:/ : /g;
		$TheLine =~ s/-/ - /g;
		$TheLine =~ s/\(\s+/ \(/g;
		$TheLine =~ s/\s+\}/\} /g;
		$TheLine =~ s/\{\s+/ \{/g;
		$TheLine =~ s/\s+\)/\) /g;
		$TheLine =~ s/=/ = /g;
		$TheLine =~ s/,/ , /g;
		$TheLine =~ s/\t/ /g;
		$TheLine =~ s/  +/ /g;
		$TheLine =~ s/^ //g;
		$TheLine =~ s/^\n+//g;
		$TheLine =~ s/ \n/\n/g;
		$TheLine =~ s/\n\n+/\n/g;
		$TheLine =~ s/\]\n/\] \n/g;

		@AddingLine = split(/ /, $TheLine);
		unless(@AddingLine == 0){
			push(@OriginalNEXUSstore, @AddingLine)
			}
		}

	#remove comments in square brackets from file
	#and make e.g.  'name one' a single entry
	@BracketNEXUSstore = @OriginalNEXUSstore;
	&ClearSquareBrackets;
	@CompoundNEXUSstore = @NewBracketNEXUSstore;
	&Compounder;
	@OriginalNEXUSstore = @CompoundNEXUSstore;

	#create a nexus store with special format for quick data searches
	@FormatNEXUSstore = ();
	@NewBracketNEXUSstore = ();
	@BracketNEXUSstore = ();

	for ($i = 0; $i <= $#MacPCNexus; $i++){
		$TheLine = $MacPCNexus[$i];
		chomp $TheLine;
		$TheLine =~ s/\\/ \\ /g;
		$TheLine =~ s/\// \/ /g;
		$TheLine =~ s/\[/ [ /g;
		$TheLine =~ s/\]/ ] /g;
		$TheLine =~ s/;/ ; /g;
		$TheLine =~ s/\s/ /g;
		$TheLine =~ s/\"/ \" /g;
		$TheLine =~ s/=/ = /g;
		$TheLine =~ s/\-/ \- /g;
		$TheLine =~ s/\(/ /g;
		$TheLine =~ s/\)/ /g;
		$TheLine =~ s/:/ /g;
		$TheLine =~ s/(dollo|irrev) *\. *up/$1\.up/ig;
		$TheLine =~ s/(dollo|irrev) *\. *dn/$1\.dn/ig;
		$TheLine =~ s/\,/ \, /g;
		$TheLine =~ s/\# +/\#/g;
		$TheLine =~ s/\* +/\*/g;
		$TheLine =~ s/  +/ /g;
		$TheLine =~ s/^ //g;

		@AddingLine = split(/\s/, $TheLine);
		unless(@AddingLine == 0){
			push(@FormatNEXUSstore, @AddingLine)
			}
		}

	#remove comments in square brackets from file
	#and make e.g.  'name one' a single entry
	@BracketNEXUSstore = @FormatNEXUSstore;
	&ClearSquareBrackets;
	@CompoundNEXUSstore = @NewBracketNEXUSstore;
	&Compounder;
	@FormatNEXUSstore = @CompoundNEXUSstore;
	}
#end of sub to get the NEXUS file input in basic form
#########################################################


#########################################################
#sub to remove comments from Nexus file
sub ClearSquareBrackets {
	@NewBracketNEXUSstore = ();
	$CommentCheck = 0;
	for($i = 0; $i <= $#BracketNEXUSstore; $i++){
		$testbit = $BracketNEXUSstore[$i];
		if($testbit =~ /\[/){$CommentCheck = 1}
		elsif($testbit =~ /\]/){$CommentCheck = 2}
		if($CommentCheck == 0){
			#add newline if needed
			if($testbit eq "\n"){
				if($NewBracketNEXUSstore[$#NewBracketNEXUSstore] !~ /\n$/){
					$NewBracketNEXUSstore[$#NewBracketNEXUSstore] = $NewBracketNEXUSstore[$#NewBracketNEXUSstore] . "\n";
					}
				}

			else{
				push(@NewBracketNEXUSstore, $testbit)
				}
			}
		if($CommentCheck == 2){$CommentCheck = 0}
		}

	$MacCladeCheck = 0;
	@Tempish = ();
	for($t = 0; $t <= $#NewBracketNEXUSstore; $t++){
		$x = $NewBracketNEXUSstore[$t];
		$x =~ s/^\n+//g;
		$y = $NewBracketNEXUSstore[$t + 1];
		$y =~ s/^\n+//g;
		$d = $x;
		$d =~ tr/A-Z/a-z/;
		$e = $y;
		$e =~ tr/A-Z/a-z/;

		if(($d =~ /^begin$/) and ($e =~ /^macclade$/)){
			$MacCladeCheck = 1;
			}
		if(($d =~ /^end$/) and ($e =~ /^;$/)){
			$MacCladeCheck = 0;
			}
		if($MacCladeCheck == 0){
			push(@Tempish, $x);
			}
		elsif($MacCladeCheck == 1){
			next;
			}
		}
	@NewBracketNEXUSstore = @Tempish;
	}
#end of sub to remove comments from Nexus file
#########################################################


#########################################################
#sub to make compound names a single entry
sub Compounder {
	@Tempy = ();
	for($SJ = 0; $SJ <= $#CompoundNEXUSstore; $SJ++){
		$TestLine = $CompoundNEXUSstore[$SJ];
		if(($TestLine =~ /^\'/) or ($TestLine =~ /^\*\'/)){
			unless($TestLine =~ /\'{1}$|\'{3}$|\'{5}$/){
				$AddingLine = $TestLine;
				$SJ++;
				$TestLine = $CompoundNEXUSstore[$SJ];
				$AddingLine = $AddingLine . " " . $TestLine;
				until($TestLine =~ /\'{1}$|\'{3}$|\'{5}$/){
					$SJ++;
					$TestLine = $CompoundNEXUSstore[$SJ];
					$AddingLine = $AddingLine . " " . $TestLine;
					}
				$TestLine = $AddingLine;
				}
			}
		push(@Tempy, $TestLine);
		}
	@CompoundNEXUSstore = @Tempy;
	}
#end of sub to make compound names a single entry
#########################################################


#########################################################
#sub to get taxon labels
sub GetTaxLabels {
	$foundTaxLabels = 0;
	$NonDefaultTaxa = 0;
	@TempTemp = ('taxset');

	#make default labels
	for($i = 1; $i <= $NTAX; $i++){
		push(@TaxLabels, $i);
		push(@TempTemp, 1);
		$TaxLabelTotalArray{"$i"} = [('taxon', $i)];
		}
	$TaxLabelTotalArray{"all"} = [@TempTemp];

	#look for taxlabels
	$Newton = 0;
	for($i = 0; $i <= $#UserNEXUSstore; $i++){
		$testbit = $UserNEXUSstore[$i];

		#check for newtaxa block
		if($testbit =~ /^newtaxa$/i){
			@TaxLabels = ();
			%TaxLabelTotalArray = ();
			$NonDefaultTaxa = 0;
			@TempTemp = ('taxset');
			for($zz = 1; $zz <= $NTAX; $zz++){
				push(@TaxLabels, $zz);
				$TaxLabelTotalArray{"$zz"} = [('taxon', $zz)];
				push(@TempTemp, 1);
				}
			$TaxLabelTotalArray{"all"} = [@TempTemp];
			$Newton = 1;
			next;
			}

		#load in taxlabels
		if($testbit =~ /^taxlabels$/i){
			$TaxCounter = 0;
			$TaxCounter2 = 1;
			$i++;
			$testbit = $UserNEXUSstore[$i];

			until($testbit =~ /\;/){
				chomp($testbit);

				unless($testbit eq "_"){
					$foundTaxLabels = 1;
					$NonDefaultTaxa = 1;
					$TaxLabels[$TaxCounter] = $testbit;
					unless(exists($TaxLabelTotalArray{"$testbit"})){
						$TaxLabelTotalArray{"$testbit"} = [('taxon', $TaxCounter2)];
						}
					}
				$TaxCounter++;
				$TaxCounter2++;
				$i++;
				$testbit = $UserNEXUSstore[$i];
				}
			if($Newton == 1){
				last;
				}
			}
		}
	}
#end of sub to get taxon labels
#########################################################


#########################################################
#sub to get Character labels
sub GetCharLabels {
	$foundCharLabels = 0;
	$NonDefaultChars = 0;

	#make default labels
	for($j = 1; $j <= $NCHAR; $j++){
		push(@CharLabels, "_");
		$TotalCharLabels{"$j"} = [('char', $j)];
		}

	#now look for charlabel block
	for($i = 0; $i <= $#UserNEXUSstore; $i++){
		$testbit = $UserNEXUSstore[$i];
		chomp($testbit);
		if($testbit =~ /^charlabels$/i){
			$foundCharLabels = 1;
			$i++;
			$L = 0;
			$testbit = $UserNEXUSstore[$i];
			last;
			}
		}

	if($foundCharLabels == 1){
		until($testbit =~ /\;/){
			chomp($testbit);
			unless($testbit =~ /_/){
				$NonDefaultChars = 1;
				$CharLabels[$L] = $testbit;
				unless($testbit eq "_"){
					$temp = $L + 1;
					$TotalCharLabels{"$testbit"} = [('char', $temp)];
					}
				}
			$i++;
			$L++;
			$testbit = $UserNEXUSstore[$i];
			}
		}
	}
#end of sub to get Charcter labels
#########################################################


#########################################################
#sub to get CharState List
sub GetCharStateLabels {
	$foundCharStateLabels = 0;

	for($i = 0; $i <= $#UserNEXUSstore; $i++){
		$testbit = $UserNEXUSstore[$i];
		chomp($testbit);
		if($testbit =~ /^charstatelabels$/i){
			$foundCharStateLabels = 1;
			$i++;
			$L = 0;
			$testbit = $UserNEXUSstore[$i];
			chomp($testbit);
			last;
			}
		}

	if($foundCharStateLabels == 1){
		$breaker = 0;
		until($testbit =~ /^\;$/){
			@templabeldata = ('char', $testbit);
			$TempName = $testbit;
			$TempNumber = $testbit;
			$startstates = 0;
			$i++;
			$testbit = $UserNEXUSstore[$i];
			chomp($testbit);
			#find data for this char...

			until($testbit =~ /^,$|^\;$/) {

				if($startstates == 0){
					if($testbit =~ /^\/$/){
						$startstates = 1;
						$i++;
						$testbit = $UserNEXUSstore[$i];
						chomp($testbit);
						}

					else{
						$TempName = $testbit;
						$i++;
						$testbit = $UserNEXUSstore[$i];
						chomp($testbit);
						next;
						}
					}

				if($startstates == 1){
					if($testbit =~ /^\;$/){$breaker = 1; last}
					if($testbit =~ /^\,$/){next}
					push(@templabeldata, $testbit);
					}
				$i++;
				$testbit = $UserNEXUSstore[$i];
				chomp($testbit);
				if($testbit =~ /^\;$/){$breaker = 1}
				}

			#add to lists:
			#as a char name, if applicable...
			unless(($TempName eq "_") or ($TempName == $TempNumber)){
				$NonDefaultChars = 1;
				$CharLabels[$TempNumber - 1] = $TempName;
				$TotalCharLabels{"$TempName"} = [@templabeldata];
				}
			#and as a char number...
			$TotalCharLabels{"$TempNumber"} = [@templabeldata];

			unless($breaker == 1){
				$i++;
				$testbit = $UserNEXUSstore[$i];
				chomp($testbit);
				}
			}
		}
	}
#end of sub to get CharState List
#########################################################


#########################################################
#sub to check the file is in NEXUS format
sub NexusCheck {
	$NexusCheck = 0;
	for($i = 0; $i <= $#UserFORMATstore; $i++){
		$testbit = $UserFORMATstore[$i];
		if($testbit =~ /\#nexus$/i){
			$NexusCheck = 1;
			last;
			}
		}
	}
#end of sub to check the file is in NEXUS format
#########################################################


#########################################################
#sub to check for datatype
sub DataTypeCheck {
	$DataType = "standard";
	$record = 0;
	@temp = ();

	for($gF = 0; $gF <= $#UserFORMATstore; $gF++){
		$testbit = $UserFORMATstore[$gF];
		if($testbit =~ /^format$/i){
			$record = 1;
			}
		if($testbit =~ /^;$/i){
			$record = 0;
			}
		if($record == 1){
			push(@temp, $testbit);
			}
		}

	if(@temp > 0){
		for($gF = 0; $gF <= $#temp; $gF++){
			$DataSearch = $temp[$gF];
			if($DataSearch =~ /^datatype$/i){
				$gF += 2;
				$DataType = $temp[$gF];
				$DataType =~ tr/A-Z/a-z/;
				last;
				}
			}
		}

	if($DataType ne "standard"){
		$RespectCase = 0;
		}

	#make settings for various data types...
	if($DataType eq "dna"){
		@tempFrom = ('R','Y','M','K','S','W','H','B','V','D','N');
		@tempTo = ('{AG}','{CT}','{AC}','{GT}','{CG}','{AT}','{ACT}','{CGT}','{ACG}','{AGT}','{ACGT}');
		push(@From, @tempFrom);
		push(@To, @tempTo);

		@ExtraSymbols = ('A','G','C','T');
		if($CharSymbolsCheck == 0){
			@SymbolsII = @ExtraSymbols;
			}
		else{push(@SymbolsII, @ExtraSymbols);}

		}

	elsif($DataType eq "rna"){
		@tempFrom = ('R','Y','M','K','S','W','H','B','V','D','N');
		@tempTo = ('{AG}','{CU}','{AC}','{GU}','{CG}','{AU}','{ACU}','{CGU}','{ACG}','{AGU}','{ACGU}');
		push(@From, @tempFrom);
		push(@To, @tempTo);

		@ExtraSymbols = ('A','G','C','U');
		if($CharSymbolsCheck == 0){
			@SymbolsII = @ExtraSymbols;
			}
		else{
			push(@SymbolsII, @ExtraSymbols);
			}
		}

	elsif($DataType eq "nucleotide"){
		@tempFrom = ('T','R','Y','M','K','S','W','H','B','V','D','N');
		@tempTo = ('U','{AG}','{CU}','{AC}','{GU}','{CG}','{AU}','{ACU}','{CGU}','{ACG}','{AGU}','{ACGU}');
		push(@From, @tempFrom);
		push(@To, @tempTo);

		@ExtraSymbols = ('A','G','C','T');
		if($CharSymbolsCheck == 0){
			@SymbolsII = @ExtraSymbols;
			}
		else{
			push(@SymbolsII, @ExtraSymbols);
			}
		}

	elsif($DataType eq "protien"){
		@tempFrom = ('B','Z');
		@tempTo = ('{DN}','{EQ}');
		push(@From, @tempFrom);
		push(@To, @tempTo);

		@ExtraSymbols = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*');
		if($CharSymbolsCheck == 0){
			@SymbolsII = @ExtraSymbols;
			}
		else{
			push(@SymbolsII, @ExtraSymbols);
			}
		}
	}

#end of sub to check for datatype
#########################################################


#########################################################
#sub to check for data layout
sub LayoutCheck {
	$Inverted = 0;
	$Interleaf = 0;
	$RespectCase = 0;

	for($i = 0; $i <= $#UserFORMATstore; $i++){
		$testbit = $UserFORMATstore[$i];
		if($testbit =~ /^interleave$/i){
			$Interleaf = 1;
			}
		if($testbit =~ /^transpose$/i){
			$Inverted = 1;
			}
		if($testbit =~ /^respectcase$/i){
			$RespectCase = 1;
			}
		}

	}
#end of sub to check for data layout
#########################################################


#########################################################
#sub to check for NTAXA
sub NumberTAXA {
	$Newton = 0;

	for($i = 0; $i <= $#UserFORMATstore; $i++){
		$testbit = $UserFORMATstore[$i];

		if($testbit =~ /^newtaxa$/i){
			$Newton = 1;
			}

		if($testbit =~ /^ntax$/i){
			$NTAX = $UserFORMATstore[$i + 2];
			if($Newton == 1){last}
			}
		}
	}

#########################################################


#########################################################
#sub to check for NCHAR
sub NumberCHARS {
	$MatrixLabels = "left";

	for($i = 0; $i <= $#UserFORMATstore; $i++){
		$testbit = $UserFORMATstore[$i];

		if($testbit =~ /^labels$/i){
			$tempton = $UserFORMATstore[$i + 2];
			if($tempton =~ /^right$/i){
				$MatrixLabels = "right"}
			if($tempton =~ /^no$/i){
				$MatrixLabels = "no"}
			}

		if($testbit =~ /^nchar$/i){
			$NCHAR = $UserFORMATstore[$i + 2];
			}
		}
	}
#end of sub to check for NCHAR
#########################################################


#########################################################
#sub to check for Missing
sub MissingSymbol {
	$Missing = "?";
	$GapMode = "missing";
	$FoundGap = 0;

	for($i = 0; $i <= $#UserNEXUSstore; $i++){
		$testbit = $UserNEXUSstore[$i];
		chomp($testbit);
		if(($testbit =~ /^missing$/i) and ($UserNEXUSstore[$i - 1] !~ /^=$/i)){
			$Missing = $UserNEXUSstore[$i + 2];
			chomp($Missing);
			}
		if($testbit =~ /^matchchar$/i){
			$Matchchar = $UserNEXUSstore[$i + 2];
			chomp($Matchchar);
			}
		if($testbit =~ /^gap$/i){
			$Gap = $UserNEXUSstore[$i + 2];
			chomp($Gap);
			$FoundGap = 1;
			}
		if($testbit =~ /^gapmode$/i){
			$GapMode = $UserNEXUSstore[$i + 2];
			$GapMode =~ tr/A-Z/a-z/;
			chomp($GapMode);
			if($GapMode eq "indel"){
				$GapMode = "missing";
				}
			}
		}
	}
#end of sub to check for Missing
#########################################################


#########################################################
#sub to check for Character Symbols
sub CharSymbols {
	$CharSymbolsCheck = 0;

	for($i = 0; $i <= $#UserNEXUSstore; $i++){
		$testbit = $UserNEXUSstore[$i];
		chomp($testbit);
		if($testbit =~ /^symbols$/i){
			$CharSymbolsCheck = 1;
			last
			}
		}

	#if symbols defined...
	if($CharSymbolsCheck == 1){
		$i += 2;
		$testbit = $UserNEXUSstore[$i];
		chomp($testbit);
		if($testbit =~ /^\"$/){
			$i ++;
			$testbit = $UserNEXUSstore[$i];
			chomp($testbit);
			}

		else{$testbit =~ s/^\"//g};

		until($testbit =~ /\"$/){
			if(length($testbit) > 1){
				@DividedBit = split(//, $testbit);
				if($DividedBit[$#DividedBit] =~ /\"/){
					pop(@DividedBit);
					push(@SymbolsII, @DividedBit);
					last;
					}
				else{
					push(@SymbolsII, @DividedBit)
					}
				}
			elsif(length($testbit) == 1){
				if($testbit eq "\""){
					last;
					}

				push(@SymbolsII, $testbit)
				}
			$i++;
			$testbit = $UserNEXUSstore[$i];
			chomp($testbit);
			}

		unless($testbit =~ /^\"$/){
			chop($testbit);
			@DividedBit = split(//, $testbit);
			push(@SymbolsII, @DividedBit)
			}
		}

	#default setting...
	if($CharSymbolsCheck == 0){
		@SymbolsII = (0,1);
		}

	}
#end of sub to check for Character Symbols
#########################################################


#########################################################
#sub to check for Equated symbols
sub EquateSymbols {
	$EquateSymbolsCheck = 0;
	@TempChopper = ();
	$DoIcheck = 0;

	for($StartCol = 0; $StartCol <= $#UserNEXUSstore; $StartCol++){
		$LineIn = $UserNEXUSstore[$StartCol];
		chomp($LineIn);
		if($LineIn =~ /^equate$/i){
			$EquateSymbolsCheck = 1;
			last
			}
		}

	if($EquateSymbolsCheck == 1){
		$StartCol += 2;
		$LineIn = $UserNEXUSstore[$StartCol];
		chomp($LineIn);

		if($LineIn =~ /^\"$/){
			$StartCol ++;
			$LineIn = $UserNEXUSstore[$StartCol];
			chomp($LineIn);
			@TempChopper = split(//, $LineIn);
			}
		else{
			$LineIn =~ s/^\"//g;
			@TempChopper = split(//, $LineIn);
			};

		while(1){
			#get from
			if($TempChopper[0] =~ /\(|\{/){
				$TempChar = shift(@TempChopper);
				&PolymorphReader;

				$testbit = $TempChar;
				$TempChar = "";
				}

			else{
				$testbit = shift(@TempChopper);
				}
			push(@From, $testbit);

			#get to
			$StartCol += 2;
			$testbit = $UserNEXUSstore[$StartCol];
			chomp($testbit);
			if($testbit =~ /\"$/){last}
			@TempChopper = split(//, $testbit);

			if($TempChopper[0] =~ /\(|\{/){
				$TempChar = shift(@TempChopper);
				&PolymorphReader;
				$testbit = $TempChar;
				$TempChar = "";
				}

			else{
				$testbit = shift(@TempChopper);
				}
			push(@To, $testbit);

			$StartCol++;
			$testbit = $UserNEXUSstore[$StartCol];
			chomp($testbit);
			if($testbit =~ /\"$/){last}
			@TempChopper = split(//, $testbit);
			}

		unless($testbit =~ /^\"$/){
			chop($testbit);
			@TempChopper = split(//, $testbit);
			if($TempChopper[0] =~ /\(|\{/){
				$TempChar = shift(@TempChopper);
				&PolymorphReader;
				$testbit = $TempChar;
				$TempChar = "";
				}

			else{
				$testbit = shift(@TempChopper);
				}

			push(@To, $testbit);
			}
		}
	}
#end of sub to check for Equated symbols
#########################################################


#########################################################
#sub to set create-date
sub CreateDate {
	($CreateSec, $CreateMin, $CreateHour, $CreateDate, $CreateMonth, $CreateYear, $CreateDay, $Used, $Used) = localtime(time);
	@Days = ('Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday');
	@Months = ('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December');
	$CreateYear = $CreateYear + 1900;
	if(($CreateDate =~ /1$/)and($CreateDate !~ /11/)){$CreateDate = $CreateDate . "st"}
	elsif(($CreateDate =~ /2$/)and($CreateDate !~ /12/)){$CreateDate = $CreateDate . "nd"}
	elsif(($CreateDate =~ /3$/)and($CreateDate !~ /13/)){$CreateDate = $CreateDate . "rd"}
	else{$CreateDate = $CreateDate . "th"}
	}
#sub to set create-date
#########################################################


#########################################################
#Sub to change the list seperator, put here because
#the $" character fucks up Textpad's formatting!
sub NewListSplit {$" = $Change} #"
#########################################################

#########################################################
#The end
