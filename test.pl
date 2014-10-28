#!/usr/bin/perl
use strict;
use warnings;

package genomeID;
use MIME::Base64;
use Bio::DB::Sam;
use List::MoreUtils 'first_index';
use Data::Dumper;

# define input variables
my $VERSION = 1;
my $file; my $type; my $noise; my $sex;
my $prefix; my $ucn; my $ref; my $baq;
my $guess_hg; my $ref_hg; my $sampleName;
my $readCol = 1;

my $Y_marker_id;
my $max_depth;
my $Y_marker_name;

# define static constants
my %bit_loc_static; my %pengelly_static; my %autoSomes_static;
my %allosomes_static; my %ref_allel_static;

my %bit_loc; my %pengelly; my %autoSomes; my %allosomes;
my %ref_allel;

# define bit holders
my $MADIB; my @AMB; my @SMB; 
my $GVB; my @XMB;

# some precompiled regexes
my $homo1 = qr/^(1[\/|]1)[^\/]?/;
my $homo2 = qr/^0[\/|]0[^\/]?/;
my $misR = qr/^(\.[\/|]\.)/;

my %chr1_lengths = ( '248956422'=>'hg38', '249250621'=>'hg19');

sub regen{
# define marker relative bit locations
%bit_loc_static = (
					#autosome bit markers
	'hg19' => ['1:45973928'=> 1, '1:67861520'=> 2, '1:158582646'=> 3, '1:179520506'=> 4,
					'1:209968684'=> 5, '1:228431095'=> 6, '2:75115108'=> 7, '2:169789016'=> 8,
					'2:215820013'=> 9, '2:227896976'=> 10, '2:49381585'=> 11, '3:4403767'=> 12,
					'3:45989044'=> 13, '3:148727133'=> 14, '4:5749904'=> 15, '4:86915848'=> 16,
					'5:13719022'=> 17, '5:40981689'=> 18, '5:55155402'=> 19, '5:82834630'=> 20,
					'5:138456815'=> 21, '5:171849471'=> 22, '7:34009946'=> 23, '7:48450157'=> 24,
					'7:100804140'=> 25, '7:151254175'=> 26, '8:94935937'=> 27, '9:27202870'=> 28,
					'9:77415284'=> 29, '9:100190780'=> 30, '10:69926097'=> 31, '10:100219314'=> 32,
					'11:6629665'=> 33, '11:16133413'=> 34, '11:30255185'=> 35, '12:993930'=> 36,
					'12:52200742'=> 37, '13:25466955'=> 38, '13:39433606'=> 39, '14:50769717'=> 40,
					'14:76045858'=> 41, '16:70303580'=> 42, '17:7192091'=> 43, '17:42449789'=> 44,
					'17:71197748'=> 45, '18:21413869'=> 46, '19:10267077'=> 47, '19:33353464'=> 48,
					'19:55441902'=> 49, '20:6100088'=> 50, '20:19970705'=> 51, '20:35865054'=> 52,
					'20:52786219'=> 53, '21:44323590'=> 54, '22:21141300'=> 55, '22:37469591'=> 56,
					'6:56471402'=> 57, '6:146755140'=> 58,

					#allosome bit markers
					'X:4066743'=>1,    'X:5616964'=>13,  'X:11882557'=>25, 'X:18151389'=>37,
					'X:25779585'=>3, 'X:34091679'=>15, 'X:35632777'=>27, 'X:40871567'=>39,
					'X:42848162'=>5, 'X:94465580'=>17, 'X:116521416'=>29,'X:139920048'=>41,
					'X:47686005'=>7, 'X:95114611'=>19, 'X:119826608'=>31,'X:144091597'=>43,
					'X:67200648'=>9, 'X:104268001'=>21,'X:120864696'=>33,'X:145069663'=>45,
					'X:78397143'=>11,'X:106598707'=>23,'X:126325138'=>35,'X:151150133'=>47,

					'Y:14902417'=>'1:4:A0-L994','Y:6788191'=>'2:2:A0-V148','Y:18255814'=>'3:4:A0-V164','Y:16292866'=>'4:2:A0-V148',
					'Y:21739754'=>'5:3:A1-M31','Y:6846482'=>'6:3:A1-M31','Y:2710154'=>'7:2:A1234-V168','Y:17947672'=>'8:2:A1234-V168',
					'Y:4898665'=>'9:2:A1234-V168','Y:21722268'=>'10:6:A2-P3','Y:18571026'=>'11:7:A2-M6','Y:21881115'=>'12:8:A2-P28',
					'Y:6845936'=>'13:5:A2-V50','Y:7589303'=>'14:3:A234-P108','Y:21722098'=>'15:7:A3-M13','Y:21925500'=>'16:6:A3-M144',
					'Y:15029492'=>'17:7:A3-M13','Y:21729839'=>'18:6:A3-M28','Y:21866840'=>'19:3:A4-M42','Y:23497067'=>'20:7:B-M192',
					'Y:21932926'=>'21:12:B-G1','Y:21774180'=>'22:13:B-M109','Y:21767959'=>'23:7:B-112','Y:21869519'=>'24:7:B-M150',
					'Y:14869076'=>'25:6:B-M182','Y:15014262'=>'26:8:B-192','Y:15437063'=>'27:10:B-218','Y:2657411'=>'28:9:B-P6',
					'Y:2711408'=>'29:9:B-M17014','Y:6982169'=>'30:9:B-M17014','Y:2844095'=>'31:10:B-M7583','Y:6753433'=>'32:10:B-M7583',
					'Y:7606120'=>'33:11:B-M8035','Y:7769194'=>'34:11:B-M8035','Y:6790204'=>'35:10:B-M8351','Y:2736732'=>'36:9:B-M8495',
					'Y:8156900'=>'37:9:B-M8495','Y:6768265'=>'38:10:B-P6','Y:14850579'=>'39:5:B-P85','Y:23120337'=>'40:8:B-V75',
					'Y:7599027'=>'41:10:B-V227','Y:2798459'=>'42:4:B-M60','Y:17201891'=>'43:11:B-V62','Y:17421349'=>'44:9:B-V75',
					'Y:6859558'=>'45:8:B-V78','Y:21866491'=>'46:9:C-M8','Y:15576203'=>'47:8:C-M208','Y:15575780'=>'48:8:C-M210',
					'Y:15437564'=>'49:6:C-M130','Y:15437333'=>'50:7:C-M217','Y:21742158'=>'51:7:C-M38','Y:21749881'=>'52:9:C-M48',
					'Y:7291534'=>'53:9:C-M8','Y:6845955'=>'54:9:C-V20','Y:14813991'=>'55:5:CDEF-M168','Y:14197867'=>'56:6:CF-P143',
					'Y:2828425'=>'57:7:D-M174','Y:24464597'=>'58:11:D-JST022457','Y:21765281'=>'59:9:D-L1375','Y:21930287'=>'60:10:D-M125',
					'Y:21872738'=>'61:8:D-M55','Y:14851526'=>'62:9:D-N1','Y:8852858'=>'63:13:E-L485','Y:14636457'=>'64:12:E-L514',
					'Y:21896261'=>'65:9:E-M33','Y:21778998'=>'66:12:E-M281','Y:22744939'=>'67:15:E-M293','Y:21740450'=>'68:9:E-M33',
					'Y:21741703'=>'69:12:E-M35','Y:21752644'=>'70:10:E-M44','Y:6822948'=>'71:14:E-M521','Y:21890177'=>'72:8:E-M75',
					'Y:21892572'=>'73:15:E-M81','Y:21778998'=>'74:7:E-M96','Y:21083420'=>'75:8:E-P147','Y:14159846'=>'76:9:E-P177',
					'Y:14060308'=>'77:10:E-P2','Y:16253694'=>'78:14:E-U175','Y:16374426'=>'79:17:E-U181','Y:17294958'=>'80:15:E-U209',
					'Y:21646058'=>'81:16:E-U290','Y:6859957'=>'82:15:E-V22','Y:6932821'=>'83:17:E-V32','Y:6818291'=>'84:11:E-V38',
					'Y:6932007'=>'85:13:E-V6','Y:17664771'=>'86:13:E-V68','Y:7729622'=>'87:14:E-Z830','Y:14197977'=>'88:17:N-P189',
					'Y:21917313'=>'89:7:F-M89','Y:6345148'=>'90:13:G-PF3146','Y:15027529'=>'91:9:G-M201','Y:22741740'=>'92:10:G-M285',
					'Y:23243942'=>'93:10:G-M285','Y:15027433'=>'94:11:G-M377','Y:22072097'=>'95:10:G-P287','Y:15058878'=>'96:16:G-U1',
					'Y:23973594'=>'97:11:G-P15','Y:23496331'=>'98:13:G-APT','Y:21764431'=>'99:11:H-M282','Y:21753199'=>'100:12:H-M52',
					'Y:14640715'=>'101:13:I-M253','Y:14847792'=>'102:12:I-M170','Y:21717307'=>'103:16:I-M223','Y:15591446'=>'104:16:I-M227',
					'Y:15022707'=>'105:13:I-M253','Y:15023364'=>'106:12:I-258','Y:21865821'=>'107:16:I-M26','Y:15426005'=>'108:17:I-P109',
					'Y:16354708'=>'109:12:I-U179','Y:21225770'=>'110:11:IJ-P126','Y:5359116'=>'111:14:J-M12','Y:5274763'=>'112:15:J-L620',
					'Y:6932171'=>'113:17:J-L817','Y:15541333'=>'114:14:J-L287','Y:9879005'=>'115:19:J-L858','Y:21926112'=>'116:14:J-M102',
					'Y:7583480'=>'117:14:J-M12','Y:21716366'=>'118:18:J-M158','Y:15587509'=>'119:15:J-M205','Y:15018459'=>'120:15:J-M241',
					'Y:22741818'=>'121:13:J-M267','Y:15467785'=>'122:18:J-M319','Y:2751678'=>'123:14:J-M410','Y:17419934'=>'124:13:J-M497',
					'Y:21904023'=>'125:19:J-M92','Y:19179335'=>'126:12:J-P209','Y:4065584'=>'127:16:J-PF5116','Y:4929023'=>'128:18:J-Z1850',
					'Y:16359006'=>'129:17:J-Z1865','Y:22797639'=>'130:19:J-Z1884','Y:14555305'=>'131:19:J-Z467','Y:21730257'=>'132:11:KLT-M9',
					'Y:14790163'=>'133:1:ZZ-123','Y:21739646'=>'134:17:L-M27','Y:2888252'=>'135:17:L-M357','Y:21757316'=>'136:17:L-M76',
					'Y:21866424'=>'137:16:M-M106','Y:17309267'=>'138:16:M-M353','Y:2744628'=>'139:16:M-M4','Y:9004921'=>'140:16:M-P336',
					'Y:19431608'=>'141:17:N-L729','Y:21741755'=>'142:19:N-M178','Y:15471925'=>'143:15:NO-M214','Y:7546726'=>'144:20:O-JST002611',
					'Y:21907097'=>'145:19:O-M110','Y:21762685'=>'146:18:O-M119','Y:21764674'=>'147:17:O-M122','Y:14924869'=>'148:20:O-M188',
					'Y:22739301'=>'149:18:O-M268','Y:2821786'=>'150:18:O-M324','Y:21868672'=>'151:19:O-M110','Y:4078217'=>'152:21:O-M7',
					'Y:21900849'=>'153:22:O-M88','Y:21938444'=>'154:21:O-M95','Y:21889536'=>'155:21:O-N6','Y:14001232'=>'156:20:O-P164',
					'Y:7568568'=>'157:16:O-P186','Y:19179463'=>'158:18:O-P197','Y:18560005'=>'159:18:O-P200','Y:14495243'=>'160:18:O-M268',
					'Y:15999244'=>'161:22:O-M117','Y:21285962'=>'162:20:O-PK4','Y:7014317'=>'163:18:Q-L472','Y:23292782'=>'164:21:Q-L54',
					'Y:15574102'=>'165:19:Q-M346','Y:21907394'=>'166:21:Q-M120','Y:21889818'=>'167:20:Q-M25','Y:21733231'=>'168:24:Q-M19',
					'Y:15018582'=>'169:16:Q-M242','Y:21866664'=>'170:20:Q-M25','Y:19096363'=>'171:23:Q-M3','Y:4925637'=>'172:18:Q-MEH2',
					'Y:21867787'=>'173:15:QR-M45','Y:4929885'=>'174:19:R-L295','Y:14641193'=>'175:24:R-L52','Y:19054889'=>'176:22:R-L757',
					'Y:21764501'=>'177:18:R-M124','Y:15030752'=>'178:20:R-M17','Y:22739367'=>'179:21:R-M269','Y:2887824'=>'180:18:R-M343',
					'Y:8533735'=>'181:21:R-M417','Y:4446430'=>'182:22:R-M520','Y:18656508'=>'183:20:R-P297','Y:18099054'=>'184:20:R-V69',
					'Y:14181107'=>'185:27:R-Z301','Y:15481372'=>'186:16:S-M230','Y:14001024'=>'187:16:S-P202','Y:6740377'=>'188:16:S-P79',
					'Y:19372808'=>'189:16:T-L131','Y:16019072'=>'190:15:T-L162','Y:6736443'=>'191:16:T-L208','Y:21893881'=>'192:14:T-M70'],

	'hg38' => ['1:45508256'=> 1, '1:67395837'=> 2, '1:158612856'=> 3, '1:179551371'=> 4,
					'1:209795339'=> 5, '1:228243394'=> 6, '2:74887981'=> 7, '2:168932506'=> 8,
					'2:214955289'=> 9, '2:227032260'=> 10, '2:49154446'=> 11, '3:4362083'=> 12,
					'3:45947552'=> 13, '3:149009346'=> 14, '4:5748177'=> 15, '4:85994695'=> 16,
					'5:13718913'=> 17, '5:40981587'=> 18, '5:55859574'=> 19, '5:83538811'=> 20,
					'5:139121126'=> 21, '5:172422467'=> 22, '7:33970334'=> 23, '7:48410560'=> 24,
					'7:101160859'=> 25, '7:151557089'=> 26, '8:93923709'=> 27, '9:27202872'=> 28,
					'9:74800368'=> 19, '9:97428498'=> 30, '10:68166340'=> 31, '10:98459557'=> 32,
					'11:6608435'=> 33, '11:16111867'=> 34, '11:30233638'=> 35, '12:884764'=> 36,
					'12:51806958'=> 37, '13:24892817'=> 38, '13:9532292'=> 39, '14:50302999'=> 40,
					'14:75579515'=> 41, '16:70269677'=> 42, '17:7288772'=> 43, '17:44372421'=> 44,
					'17:73201609'=> 45, '18:23833905'=> 46, '19:10156401'=> 47, '19:32862558'=> 48,
					'19:54930534'=> 49, '20:6119441'=> 50, '20:19990061'=> 51, '20:37236651'=> 52,
					'20:54169680'=> 53, '21:42903480'=> 54, '22:20787012'=> 55, '22:37073551'=> 56,
					'6:56606604'=> 57, '6:146434004'=> 58,

					#allosome bit markers
					'X:4148702'=>1,'X:5698923'=>13,'X:11864438'=>25,'X:18133269'=>37,
					'X:25761468'=>3,'X:34073562'=>15,'X:35614660'=>27,'X:41012314'=>39,
					'X:42988913'=>5,'X:95210581'=>17,'X:117387453'=>29,'X:140837883'=>41,
					'X:47826606'=>7,'X:95859612'=>19,'X:120692753'=>31,'X:145010077'=>43,
					'X:67980806'=>9,'X:105023319'=>21,'X:121730843'=>33,'X:145988145'=>45,
					'X:79141646'=>11,'X:107355477'=>23,'X:127191155'=>35,'X:151981661'=>47]);
 
# declar pengelly marker locations
%pengelly_static = ( 
	'hg19' => ['1:45973928'=> 0, '1:67861520'=> 1, '1:158582646'=> 0, '1:179520506'=> 1,
					'1:209968684'=> 0, '1:228431095'=> 0, '2:75115108'=> 0, '2:169789016'=> 1,
					'2:215820013'=> 0, '2:227896976'=> 1, '2:49381585'=> 0, '3:4403767'=> 1,
					'3:45989044'=> 0, '3:148727133'=> 0, '4:5749904'=> 1, '4:86915848'=> 0,
					'5:13719022'=> 0, '5:40981689'=> 0, '5:55155402'=> 0, '5:82834630'=> 1,
					'5:138456815'=> 0, '5:171849471'=> 0, '7:34009946'=> 0, '7:48450157'=> 1,
					'7:100804140'=> 0, '7:151254175'=> 0, '8:94935937'=> 1, '9:27202870'=> 0,
					'9:77415284'=> 0, '9:100190780'=> 1, '10:69926097'=> 0, '10:100219314'=> 1,
					'11:6629665'=> 0, '11:16133413'=> 1, '11:30255185'=> 0, '12:993930'=> 1,
					'12:52200742'=> 0, '13:25466955'=> 0, '13:39433606'=> 1, '14:50769717'=> 1,
					'14:76045858'=> 0, '16:70303580'=> 1, '17:7192091'=> 0, '17:42449789'=> 0,
					'17:71197748'=> 1, '18:21413869'=> 1, '19:10267077'=> 1, '19:33353464'=> 0,
					'19:55441902'=> 0, '20:6100088'=> 1, '20:19970705'=> 0, '20:35865054'=> 0,
					'20:52786219'=> 0, '21:44323590'=> 1, '22:21141300'=> 1, '22:37469591'=> 0,
					'6:56471402'=> 0, '6:146755140'=> 1], 

	'hg38' => ['1:45508256'=> 0, '1:67395837'=> 1, '1:158612856'=> 0, '1:179551371'=> 1,
					'1:209795339'=> 0, '1:228243394'=> 0, '2:74887981'=> 0, '2:168932506'=> 1,
					'2:214955289'=> 0, '2:227032260'=> 1, '2:49154446'=> 0, '3:4362083'=> 1,
					'3:45947552'=> 0, '3:149009346'=> 0, '4:5748177'=> 1, '4:85994695'=> 0,
					'5:13718913'=> 0, '5:40981587'=> 0, '5:55859574'=> 0, '5:83538811'=> 1,
					'5:139121126'=> 0, '5:172422467'=> 0, '7:33970334'=> 0, '7:48410560'=> 1,
					'7:101160859'=> 0, '7:151557089'=> 0, '8:93923709'=> 1, '9:27202872'=> 0,
					'9:74800368'=> 0, '9:97428498'=> 1, '10:68166340'=> 0, '10:98459557'=> 1,
					'11:6608435'=> 0, '11:16111867'=> 1, '11:30233638'=> 0, '12:884764'=> 1,
					'12:51806958'=> 0, '13:24892817'=> 0, '13:9532292'=> 1, '14:50302999'=> 1,
					'14:75579515'=> 0, '16:70269677'=> 1, '17:7288772'=> 0, '17:44372421'=> 0,
					'17:73201609'=> 1, '18:23833905'=> 1, '19:10156401'=> 1, '19:32862558'=> 0,
					'19:54930534'=> 0, '20:6119441'=> 1, '20:19990061'=> 0, '20:37236651'=> 0,
					'20:54169680'=> 0, '21:42903480'=> 1, '22:20787012'=> 1, '22:37073551'=> 0,
					'6:56606604'=> 0, '6:146434004'=> 1]);

# declare autosomal markers
%autoSomes_static = (
	'hg19' => ['1:45973928'=> 'A:G', '1:67861520'=> 'A:C', '1:158582646'=> 'C:T', '1:179520506'=> 'A:G',
					'1:209968684'=> 'A:C', '1:228431095'=> 'C:T', '2:75115108'=> 'A:G', '2:169789016'=> 'C:T',
					'2:215820013'=> 'A:G', '2:227896976'=> 'C:T', '2:49381585'=> 'A:G', '3:4403767'=> 'C:T',
					'3:45989044'=> 'G:T', '3:148727133'=> 'A:G', '4:5749904'=> 'T:C', '4:86915848'=> 'C:T',
					'5:13719022'=> 'A:C', '5:40981689'=> 'A:C', '5:55155402'=> 'C:T', '5:82834630'=> 'C:T',
					'5:138456815'=> 'C:T', '5:171849471'=> 'A:G', '7:34009946'=> 'C:T', '7:48450157'=> 'C:T',
					'7:100804140'=> 'C:T', '7:151254175'=> 'C:T', '8:94935937'=> 'C:T', '9:27202870'=> 'C:T',
					'9:77415284'=> 'A:C', '9:100190780'=> 'C:T', '10:69926097'=> 'A:G', '10:100219314'=> 'A:G',
					'11:6629665'=> 'C:T', '11:16133413'=> 'A:G', '11:30255185'=> 'C:T', '12:993930'=> 'C:T',
					'12:52200742'=> 'A:C', '13:25466955'=> 'C:T', '13:39433606'=> 'A:G', '14:50769717'=> 'A:G',
					'14:76045858'=> 'A:G', '16:70303580'=> 'C:T', '17:7192091'=> 'C:T', '17:42449789'=> 'C:T',
					'17:71197748'=> 'A:G', '18:21413869'=> 'C:T', '19:10267077'=> 'A:G', '19:33353464'=> 'A:G',
					'19:55441902'=> 'C:T', '20:6100088'=> 'C:T', '20:19970705'=> 'A:G', '20:35865054'=> 'C:T',
					'20:52786219'=> 'A:G', '21:44323590'=> 'G:T', '22:21141300'=> 'C:T', '22:37469591'=> 'A:G',
					'6:56471402'=> 'A:G', '6:146755140'=> 'A:G'], 

	'hg38' => ['1:45508256'=> 'A:G', '1:67395837'=> 'A:C', '1:158612856'=> 'C:T', '1:179551371'=> 'A:G',
					'1:209795339'=> 'A:C', '1:228243394'=> 'C:T', '2:74887981'=> 'A:G', '2:168932506'=> 'C:T',
					'2:214955289'=> 'A:G', '2:227032260'=> 'C:T', '2:49154446'=> 'A:G', '3:4362083'=> 'C:T',
					'3:45947552'=> 'G:T', '3:149009346'=> 'A:G', '4:5748177'=> 'T:C', '4:85994695'=> 'C:T',
					'5:13718913'=> 'A:C', '5:40981587'=> 'A:C', '5:55859574'=> 'C:T', '5:83538811'=> 'C:T',
					'5:139121126'=> 'C:T', '5:172422467'=> 'A:G', '7:33970334'=> 'C:T', '7:48410560'=> 'C:T',
					'7:101160859'=> 'C:T', '7:151557089'=> 'C:T', '8:93923709'=> 'C:T', '9:27202872'=> 'C:T',
					'9:74800368'=> 'A:C', '9:97428498'=> 'C:T', '10:68166340'=> 'A:G', '10:98459557'=> 'A:G',
					'11:6608435'=> 'C:T', '11:16111867'=> 'A:G', '11:30233638'=> 'C:T', '12:884764'=> 'C:T',
					'12:51806958'=> 'A:C', '13:24892817'=> 'C:T', '13:9532292'=> 'A:G', '14:50302999'=> 'A:G',
					'14:75579515'=> 'A:G', '16:70269677'=> 'C:T', '17:7288772'=> 'C:T', '17:44372421'=> 'C:T',
					'17:73201609'=> 'A:G', '18:23833905'=> 'C:T', '19:10156401'=> 'A:G', '19:32862558'=> 'A:G',
					'19:54930534'=> 'C:T', '20:6119441'=> 'C:T', '20:19990061'=> 'A:G', '20:37236651'=> 'C:T',
					'20:54169680'=> 'A:G', '21:42903480'=> 'G:T', '22:20787012'=> 'C:T', '22:37073551'=> 'A:G',
					'6:56606604'=> 'A:G', '6:146434004'=> 'A:G']);

# declare sex chromosome markers
%allosomes_static = ( 
	'hg19' => ['X:4066743'=>'G:A','X:5616964'=>'G:A','X:11882557'=>'T:C','X:18151389'=>'G:A',
				'X:25779585'=>'T:C','X:34091679'=>'T:C','X:35632777'=>'C:T','X:40871567'=>'G:A',
				'X:42848162'=>'C:T','X:94465580'=>'T:G','X:116521416'=>'G:T','X:139920048'=>'T:C',
				'X:47686005'=>'C:T','X:95114611'=>'T:C','X:119826608'=>'T:G','X:144091597'=>'G:C',
				'X:67200648'=>'T:C','X:104268001'=>'T:C','X:120864696'=>'C:A','X:145069663'=>'G:A',
				'X:78397143'=>'T:C','X:106598707'=>'A:G','X:126325138'=>'G:T','X:151150133'=>'G:T',

				'Y:14902417'=>'C:T','Y:6788191'=>'G:A','Y:18255814'=>'C:T','Y:16292866'=>'G:A',
				'Y:21739754'=>'G:C','Y:6846482'=>'G:A','Y:2710154'=>'G:A','Y:17947672'=>'G:A',
				'Y:4898665'=>'C:G','Y:21722268'=>'T:C','Y:18571026'=>'T:C','Y:21881115'=>'C:T',
				'Y:6845936'=>'T:C','Y:7589303'=>'T:C','Y:21722098'=>'G:C','Y:21925500'=>'T:C',
				'Y:15029492'=>'T:G','Y:21729839'=>'T:G','Y:21866840'=>'A:T','Y:23497067'=>'C:G',
				'Y:21932926'=>'A:G','Y:21774180'=>'C:T','Y:21767959'=>'G:A','Y:21869519'=>'C:T',
				'Y:14869076'=>'C:T','Y:15014262'=>'C:T','Y:15437063'=>'C:T','Y:2657411'=>'G:A',
				'Y:2711408'=>'C:G','Y:6982169'=>'A:C','Y:2844095'=>'G:T','Y:6753433'=>'A:G',
				'Y:7606120'=>'T:C','Y:7769194'=>'C:T','Y:6790204'=>'C:G','Y:2736732'=>'T:G',
				'Y:8156900'=>'T:G','Y:6768265'=>'G:C','Y:14850579'=>'C:T','Y:23120337'=>'G:A',
				'Y:7599027'=>'A:T','Y:2798459'=>'C:T','Y:17201891'=>'G:A','Y:17421349'=>'G:A',
				'Y:6859558'=>'C:G','Y:21866491'=>'C:T','Y:15576203'=>'C:T','Y:15575780'=>'A:T',
				'Y:15437564'=>'C:T','Y:15437333'=>'A:C','Y:21742158'=>'T:G','Y:21749881'=>'A:G',
				'Y:7291534'=>'G:T','Y:6845955'=>'G:A','Y:14813991'=>'C:T','Y:14197867'=>'G:A',
				'Y:2828425'=>'A:G','Y:24464597'=>'G:C','Y:21765281'=>'T:A','Y:21930287'=>'T:C',
				'Y:21872738'=>'T:C','Y:14851526'=>'C:G','Y:8852858'=>'C:T','Y:14636457'=>'C:T',
				'Y:21896261'=>'G:T','Y:21778998'=>'C:T','Y:22744939'=>'T:G','Y:21740450'=>'A:C',
				'Y:21741703'=>'G:C','Y:21752644'=>'G:C','Y:6822948'=>'C:T','Y:21890177'=>'G:A',
				'Y:21892572'=>'C:T','Y:21778998'=>'C:G','Y:21083420'=>'T:A','Y:14159846'=>'C:T',
				'Y:14060308'=>'A:C','Y:16253694'=>'G:A','Y:16374426'=>'C:T','Y:17294958'=>'C:T',
				'Y:21646058'=>'T:A','Y:6859957'=>'T:C','Y:6932821'=>'G:C','Y:6818291'=>'C:T',
				'Y:6932007'=>'G:C','Y:17664771'=>'A:C','Y:7729622'=>'G:A','Y:14197977'=>'G:A',
				'Y:21917313'=>'C:T','Y:6345148'=>'C:G','Y:15027529'=>'G:T','Y:22741740'=>'G:C',
				'Y:23243942'=>'C:T','Y:15027433'=>'A:G','Y:22072097'=>'G:T','Y:15058878'=>'A:G',
				'Y:23973594'=>'T:G','Y:23496331'=>'C:T','Y:21764431'=>'A:C','Y:21753199'=>'A:C',
				'Y:14640715'=>'A:G','Y:14847792'=>'A:C','Y:21717307'=>'G:A','Y:15591446'=>'C:G',
				'Y:15022707'=>'C:T','Y:15023364'=>'T:C','Y:21865821'=>'G:A','Y:15426005'=>'C:T',
				'Y:16354708'=>'G:A','Y:21225770'=>'C:G','Y:5359116'=>'C:T','Y:5274763'=>'G:A',
				'Y:6932171'=>'G:C','Y:15541333'=>'C:T','Y:9879005'=>'C:T','Y:21926112'=>'G:C',
				'Y:7583480'=>'G:T','Y:21716366'=>'G:A','Y:15587509'=>'T:A','Y:15018459'=>'G:A',
				'Y:22741818'=>'T:G','Y:15467785'=>'T:A','Y:2751678'=>'A:G','Y:17419934'=>'T:C',
				'Y:21904023'=>'T:C','Y:19179335'=>'T:C','Y:4065584'=>'G:A','Y:4929023'=>'G:A',
				'Y:16359006'=>'G:T','Y:22797639'=>'C:T','Y:14555305'=>'C:T','Y:21730257'=>'C:G',
				'Y:14790163'=>'A:G','Y:21739646'=>'C:G','Y:2888252'=>'C:A','Y:21757316'=>'T:G',
				'Y:21866424'=>'A:G','Y:17309267'=>'G:A','Y:2744628'=>'A:G','Y:9004921'=>'G:A',
				'Y:19431608'=>'A:C','Y:21741755'=>'C:T','Y:15471925'=>'T:C','Y:7546726'=>'C:T',
				'Y:21907097'=>'T:C','Y:21762685'=>'A:C','Y:21764674'=>'T:C','Y:14924869'=>'C:T',
				'Y:22739301'=>'A:G','Y:2821786'=>'G:C','Y:21868672'=>'T:C','Y:4078217'=>'C:G',
				'Y:21900849'=>'A:G','Y:21938444'=>'C:T','Y:21889536'=>'G:C','Y:14001232'=>'T:C',
				'Y:7568568'=>'C:A','Y:19179463'=>'G:T','Y:18560005'=>'T:G','Y:14495243'=>'T:C',
				'Y:15999244'=>'G:A','Y:21285962'=>'A:T','Y:7014317'=>'G:C','Y:23292782'=>'G:A',
				'Y:15574102'=>'G:A','Y:21907394'=>'T:C','Y:21889818'=>'G:T','Y:21733231'=>'T:A',
				'Y:15018582'=>'C:T','Y:21866664'=>'G:C','Y:19096363'=>'G:A','Y:4925637'=>'G:T',
				'Y:21867787'=>'G:A','Y:4929885'=>'G:A','Y:14641193'=>'C:T','Y:19054889'=>'C:T',
				'Y:21764501'=>'C:T','Y:15030752'=>'C:T','Y:22739367'=>'T:C','Y:2887824'=>'C:A',
				'Y:8533735'=>'G:A','Y:4446430'=>'T:A','Y:18656508'=>'G:C','Y:18099054'=>'C:T',
				'Y:14181107'=>'C:T','Y:15481372'=>'T:A','Y:14001024'=>'T:A','Y:6740377'=>'T:C',
				'Y:19372808'=>'C:T','Y:16019072'=>'G:C','Y:6736443'=>'C:T','Y:21893881'=>'A:C'],

	'hg38' => ['X:4148702'=>'G:A','X:5698923'=>'G:A','X:11864438'=>'T:C','X:18133269'=>'G:A',
				'X:25761468'=>'T:C','X:34073562'=>'T:C','X:35614660'=>'C:T','X:41012314'=>'G:A',
				'X:42988913'=>'C:T','X:95210581'=>'T:G','X:117387453'=>'G:T','X:140837883'=>'T:C',
				'X:47826606'=>'C:T','X:95859612'=>'T:C','X:120692753'=>'T:G','X:145010077'=>'G:C',
				'X:67980806'=>'T:C','X:105023319'=>'T:C','X:121730843'=>'C:A','X:145988145'=>'G:A',
				'X:79141646'=>'T:C','X:107355477'=>'A:G','X:127191155'=>'G:T','X:151981661'=>'G:T']);


# declare reference sex chromosome markers
%ref_allel_static = ( 
	'hg19' => ['X:4066743'=>'G','X:5616964'=>'A','X:11882557'=>'C','X:18151389'=>'A',
				'X:25779585'=>'T','X:34091679'=>'T','X:35632777'=>'C','X:40871567'=>'G',
				'X:42848162'=>'C','X:94465580'=>'T','X:116521416'=>'G','X:139920048'=>'C',
				'X:47686005'=>'T','X:95114611'=>'C','X:119826608'=>'G','X:144091597'=>'G',
				'X:67200648'=>'T','X:104268001'=>'C','X:120864696'=>'A','X:145069663'=>'G',
				'X:78397143'=>'C','X:106598707'=>'A','X:126325138'=>'G','X:151150133'=>'G',

				# hg19 ref allele at the hg38 allosome positions
				'X:4148702'=>'G','X:5698923'=>'A','X:11864438'=>'A','X:18133269'=>'C',
				'X:25761468'=>'A','X:34073562'=>'A','X:35614660'=>'G','X:41012314'=>'G',
				'X:42988913'=>'A','X:95210581'=>'C','X:117387453'=>'C','X:140837883'=>'G',
				'X:47826606'=>'','X:95859612'=>'C','X:120692753'=>'G','X:145010077'=>'T',
				'X:67980806'=>'T','X:105023319'=>'C','X:121730843'=>'T','X:145988145'=>'G',
				'X:79141646'=>'C','X:107355477'=>'T','X:127191155'=>'T','X:151981661'=>'G',

				# hg19 ref allele at the hg38 autosome positions
				'1:45508256'=> 'A', '1:67395837'=> 'C', '1:158612856'=> 'G', '1:179551371'=> 'C',
				'1:209795339'=> 'A', '1:228243394'=> 'A', '2:74887981'=> 'G', '2:168932506'=> 'T',
				'2:214955289'=> 'T', '2:227032260'=> 'T', '2:49154446'=> 'G', '3:4362083'=> 'A',
				'3:45947552'=> 'A', '3:149009346'=> 'A', '4:5748177'=> 'T', '4:85994695'=> 'G',
				'5:13718913'=> 'T', '5:40981587'=> 'T', '5:55859574'=> 'A', '5:83538811'=> 'T',
				'5:139121126'=> 'A', '5:172422467'=> 'A', '7:33970334'=> 'C', '7:48410560'=> 'T',
				'7:101160859'=> 'A', '7:151557089'=> 'G', '8:93923709'=> 'G', '9:27202872'=> 'C',
				'9:74800368'=> 'G', '9:97428498'=> 'A', '10:68166340'=> 'C', '10:98459557'=> 'T',
				'11:6608435'=> 'G', '11:16111867'=> 'C', '11:30233638'=> 'G', '12:884764'=> 'A',
				'12:51806958'=> 'T', '13:24892817'=> 'G', '13:9532292'=> 'N', '14:50302999'=> 'A',
				'14:75579515'=> 'C', '16:70269677'=> 'A', '17:7288772'=> 'A', '17:44372421'=> 'T',
				'17:73201609'=> 'T', '18:23833905'=> 'T', '19:10156401'=> 'T', '19:32862558'=> 'T',
				'19:54930534'=> 'T', '20:6119441'=> 'G', '20:19990061'=> 'G', '20:37236651'=> 'A',
				'20:54169680'=> 'A', '21:42903480'=> 'G', '22:20787012'=> 'G', '22:37073551'=> 'A',
				'6:56606604'=> 'A', '6:146434004'=> 'A'],

	'hg38' => ['X:4148702'=>'G','X:5698923'=>'A','X:11864438'=>'C','X:18133269'=>'A',
				'X:25761468'=>'T','X:34073562'=>'T','X:35614660'=>'C','X:41012314'=>'G',
				'X:42988913'=>'C','X:95210581'=>'T','X:117387453'=>'G','X:140837883'=>'C',
				'X:47826606'=>'T','X:95859612'=>'C','X:120692753'=>'G','X:145010077'=>'G',
				'X:67980806'=>'T','X:105023319'=>'C','X:121730843'=>'A','X:145988145'=>'G',
				'X:79141646'=>'C','X:107355477'=>'A','X:127191155'=>'G','X:151981661'=>'G',

				#hg38 autosomal ref alleles
				'1:45508256'=> 'G', '1:67395837'=> 'C', '1:158612856'=> 'T', '1:179551371'=> 'G',
				'1:209795339'=> 'C', '1:228243394'=> 'A', '2:74887981'=> 'A', '2:168932506'=> 'T',
				'2:214955289'=> 'G', '2:227032260'=> 'C', '2:49154446'=> 'C', '3:4362083'=> 'A',
				'3:45947552'=> 'T', '3:149009346'=> 'G', '4:5748177'=> 'T', '4:85994695'=> 'T',
				'5:13718913'=> 'T', '5:40981587'=> 'C', '5:55859574'=> 'C', '5:83538811'=> 'T',
				'5:139121126'=> 'T', '5:172422467'=> 'G', '7:33970334'=> 'C', '7:48410560'=> 'T',
				'7:101160859'=> 'C', '7:151557089'=> 'T', '8:93923709'=> 'T', '9:27202872'=> 'A',
				'9:74800368'=> 'A', '9:97428498'=> 'A', '10:68166340'=> 'T', '10:98459557'=> 'G',
				'11:6608435'=> 'C', '11:16111867'=> 'A', '11:30233638'=> 'C', '12:884764'=> 'C',
				'12:51806958'=> 'C', '13:24892817'=> 'T', '13:9532292'=> 'N', '14:50302999'=> 'G',
				'14:75579515'=> 'G', '16:70269677'=> 'G', '17:7288772'=> 'C', '17:44372421'=> 'G',
				'17:73201609'=> 'G', '18:23833905'=> 'T', '19:10156401'=> 'T', '19:32862558'=> 'G',
				'19:54930534'=> 'T', '20:6119441'=> 'A', '20:19990061'=> 'C', '20:37236651'=> 'C',
				'20:54169680'=> 'G', '21:42903480'=> 'T', '22:20787012'=> 'T', '22:37073551'=> 'G',
				'6:56606604'=> 'G', '6:146434004'=> 'G']);

$prefix = "";
$guess_hg  =0;
$max_depth = 0;
$Y_marker_id = 0;
$Y_marker_name = "";
$readCol = 9;

# define bit holders
$MADIB = ""; @AMB = ("0") x 58; @SMB = ("0") x 2;
$GVB = "";  @XMB = ("0") x 48;
}


sub generate_id{
	# requires input as hash
	my (%input) = @_;

	regen();

	# determine hg version
	if(exists $input{'hg'}){
		
		# determine if hg version is supported
		if(! exists $bit_loc_static{$input{'hg'}}){
			die "Version: $input{'hg'} is not supported\n";
		}
		$ref_hg = $input{'hg'};
		(%pengelly) = @{$pengelly_static{$input{'hg'}}}; 
		(%bit_loc) = @{$bit_loc_static{$input{'hg'}}};
		(%autoSomes) = @{$autoSomes_static{$input{'hg'}}};
		(%allosomes) = @{$allosomes_static{$input{'hg'}}};
		(%ref_allel) = @{$ref_allel_static{$input{'hg'}}};
	}

	# otherwise guess it by looking at marker alignment
	else{
		$guess_hg = 1;

		# start with a guess of hg19
		my $hg = 'hg38'; $ref_hg = $hg;
		(%pengelly) = @{$pengelly_static{ $hg }};
		(%bit_loc) = @{$bit_loc_static{ $hg }};
		(%autoSomes) = @{$autoSomes_static{ $hg }};
		(%allosomes) = @{$allosomes_static{ $hg }};
		(%ref_allel) = @{$ref_allel_static{ $hg }};
	}

	# store input flags
	$baq = (exists $input{'baq'})? $input{'baq'} : 30; 
	$noise = (exists $input{'noise'})? $input{'noise'} : 0.05;
	$sex = $input{'sex'}; 
	$file = $input{'file'}; 
	$ucn = (exists $input{'ucn'})? $input{'ucn'} : 0;
	$ref = (exists $input{'ref'})? $input{'ref'} : 0;
	$sampleName = (exists $input{'sampleName'})? $input{'sampleName'} : "";

	# determine if a sampleName was specified
	if( exists $input{'sampleName'} && $input{'type'} ne 'bam'){
		detReadCol($input{'type'});
		return $sampleName;
	}

	# for the XMB, we assume everything is
	# undefined, later we will fix this 
	# if necessary
	for(my $i=0;$i< scalar @XMB; $i +=2){
		$XMB[$i] = 1
	}
		
	# determine if BAM or VCF
	if ($input{'type'} eq 'bam'){
		# open file using sam tools
		$file = Bio::DB::Sam->new(-bam =>$file);
			 
		# determine if prefix is required
		my @chrs = $file->seq_ids;
		if($chrs[0] =~ /^([a-z]+)\d+/){
			$prefix = $1;
		}
		
		# determine hg version if not specified
		if($guess_hg){
			my $chrLength = $file->length("$prefix"."1");
			
			if(! exists $chr1_lengths{$chrLength}){
				die "Could not determine human genome version, please specify\n";
			}
		
			$ref_hg = $chr1_lengths{ $chrLength };
			(%pengelly) = @{$pengelly_static{ $ref_hg }};
			(%bit_loc) = @{$bit_loc_static{ $ref_hg }};
			(%autoSomes) = @{$autoSomes_static{ $ref_hg }};
			(%allosomes) = @{$allosomes_static{ $ref_hg }};
			(%ref_allel) = @{$ref_allel_static{ $ref_hg }};
			$guess_hg =0;
		}

		bam();
	}
	elsif($input{'type'} eq 'vcf'){
		vcf();
	}
	elsif($input{'type'} eq 'tbi'){
		tbi();
	}
	else{
		die "No supported file type was specified\n";
	}

	# generate version block
	$GVB = 0 x (6-length(sprintf("%b",$VERSION))) . sprintf("%b",$VERSION);
	
	# generate the Y marker bits if sex, and Y is present
	if($SMB[0] && $sex){
		$Y_marker_id = 0 x (9-length(sprintf("%b",$Y_marker_id) ) ) . sprintf("%b",$Y_marker_id);
		
		# change bits of the XMB
		for( my $i = 0; $i < 18; $i+=2){
			$XMB[$i] = substr($Y_marker_id,$i/2,1);
		}
		if($Y_marker_name eq 1){
			$XMB[18] = 1;
		}
		else{
			$XMB[18] = 0;
		}
	}

	# generate SIB
	# if ucn or non informative
	# sex marker bit, then this is 1
	$SMB[1] = ($ucn || !$SMB[0] )? 1: "0";
			
	my $XMB = ($sex)? join("",@XMB) : "";
	my $AMB = join("",@AMB);
	my $SMB = join("",@SMB);

	my $genString = "$MADIB$AMB$SMB$GVB$XMB";
	return encode_base64 pack 'B*', $genString;
}

=comment

Tabix Parser

=cut
sub tbi{
	my $missing_markers =0; my $only_peng = 1; my $bit_mis = 0;
	
	$prefix = "chr";

	# loop through autoSomal markers
	foreach my $key(keys %autoSomes){
		# grab marker location to query vcf file 
		my ($chr, $pos) = split(":",$key);
		my $isPengelly = $pengelly{$key};
		my $bit  = $bit_loc{$key} -1;
		my $mis = 0; my $zyg = 0;

		# query tabix file
		# need to figure out optimal way to determine
		# if prefix is used
		my $data = `tabix $file -b 2 -e 2 $chr:$pos-$pos`;
		if($data eq ""){
			$data = `tabix $file -b 2 -e 2 $prefix$chr:$pos-$pos`;
		}

		my @rows = split(/\n/,$data);
		my @col = split(/\t/,$data);
	  
		# sometimes query may give more 
		# than one row (i dunno why)
		if(scalar @rows > 1){
			foreach my $row (@rows){
				@col = split(/\t/,$row);
				if ($col[1] eq $pos){
					$data = $row;
					last;
				}
			}
		}
		
		# determine zygosity
		if(! $data eq ""){
			$zyg = ($col[$readCol] =~ $homo1 or $col[$readCol] =~ $homo2)? 0:1;
			$mis = ($col[$readCol] =~ $misR)? 1:0;
		}

		if($data ne "" && $guess_hg){
			my $continue = guessHG('key'=>$key,'ref'=>$col[3],'bam'=>0);
	 
	 		# determine if a change was detected
			if(!$continue){
				tbi();
				return;
			}
		}

		# no read was found
		if ($data eq "" || $mis){
			if($ref && !$mis){
				$AMB[$bit] = 0; next;
			}

			$AMB[$bit] = 0;
			$missing_markers++;
			$bit_mis = $bit;

			if($isPengelly){$only_peng =0;}
		}
		
		# read was found
		else{
			# append amb
			$AMB[$bit] = $zyg;
	 
			# if this is a non pengelly
			# then not only pengelly's are found
			if(!$isPengelly){$only_peng =0;}
		}
	}

	# generate MADIB
	genMADIB('miss_count'=>$missing_markers,'oPeng'=>$only_peng, 'bMis'=>$bit_mis);
	
	if(!$sex){
		return;
	}

	# loop through allosomal markers
	foreach my $key (keys %allosomes){
		# grab marker location to query vcf file
		my ($chr,$pos) = split(':', $key);
		my ($anc,$alt) = split(":", $allosomes{$key});

		# query tabix file
		my $data = `tabix $file -b 2 -e 2 $chr:$pos-$pos`;
		if($data eq ""){
			$data = `tabix $file -b 2 -e 2 $prefix$chr:$pos-$pos`;
		}
		my @rows = split(/\n/,$data);
		
		# determine if there was results,
		# if not, assign to reference if flag
		# was specified
		if($data eq ""){
			next if $chr eq 'Y';
			
			my $bit =  $bit_loc{$key} -1;
			if($ref){
				my $ref_al = $ref_allel{$key};

				# homozygous ancestral
				if($ref_al eq $anc && !$SMB[0]){
					$XMB[$bit] = 0;
				}

				# homozygous alternate
				else{
					$XMB[$bit+1] = 1;
				}
			}
			next;
		}

		# ensure that only one row was returned in from
		# the tabix query
		if(scalar @rows > 1){
			foreach my $row (@rows){
				my @col = split(/\t/,$row);
				
				if($col[1] eq $pos){
					$data = $row;
					last;
				}
			}
		}

		# build the xmb bits
		xmb_builder( 'key' => $key, 'data' => $data );
	}
}

=comment

BAM Parser

=cut
sub bam{
	my $missing_markers = 0; my $only_peng = 1; my $bit_mis =0;

	# loop through autoSomal markers
	foreach my $key (keys %autoSomes){  
		# grab marker location to query bam file
		my ($chr,$pos) = split(":",$key);
		my ($fir,$las) = split(":",$autoSomes{$key});
		my $isPengelly = $pengelly{$key};   
		my $bit = $bit_loc{$key} -1;
			  
		# get zygosity 
		my $zyg = bam_zygosity('chr'=>"$prefix$chr",'loc'=>$pos,'al1'=>$fir,'al2'=>$las);
		
		# if no reads
		if($zyg == -1){
			# append to AMB
			$AMB[$bit] = "0"; 
			$missing_markers++;
			$bit_mis = $bit;

			# if this is a missing pengelly, 
			# then not only pengelly's are read
			if($isPengelly){ $only_peng =0;}
		}
		
		# otherwise a read was found
		else{
			$AMB[$bit] = $zyg;

			# if this is a read non pengelly,
			# then not only pengelly's are read
			if(!$isPengelly){$only_peng = 0;}
		} 
	}
	genMADIB('miss_count'=>$missing_markers,'oPeng'=>$only_peng,'bMis'=>$bit_mis);

	$SMB[0] = (grep $_ eq  "$prefix"."Y", $file->seq_ids);
	
	# generate additional sex markers
	if($sex){
			# loop through allosomes
			foreach my $key (keys %allosomes){
			# grab marker information used to query bam file
			my ($chr,$pos) = split(":",$key);
			my ($anc,$alt) = split(":",$allosomes{$key});

			# determine zygosity
			my $zyg = bam_zygosity('chr'=>"$prefix$chr",'loc'=>$pos,'anc'=>$anc,'alt'=>$alt,'smb'=>$SMB[0]);

			# modify both bits
			if($chr eq 'X'){
				my $bit = $bit_loc{$key} -1;
				my ($bit1,$bit2) = split(":",$zyg);

				if(!$SMB[0]){
					$XMB[$bit] = $bit1;
				}
				$XMB[$bit+1] = $bit2;
			}
			
			# Y chromosome modify left bit
			else{
				next unless $zyg eq 1;
				my ($bit,$depth,$mark_id) = split(":",$bit_loc{$key});
				
				# ensure not from a difference branch
				if($Y_marker_name ne 1){
					my ($tree, $branch) = split("-", $mark_id);

					# check to see if tree is the same
					if( index($Y_marker_name, $tree) != -1 || $Y_marker_name eq ""  ){
						$Y_marker_name = $tree;
					}
					else{
						$Y_marker_name = 1;
					}
				}

				# we take the deepest branch 
				next unless $depth > $max_depth;
				
				$max_depth = $depth;
				$Y_marker_id = $bit;
			}

		}
	}
}


=comment
	
VCF Parser

=cut
sub vcf{
	# define useful constants
	my $missing_markers = 0;
	my $only_peng  = 1; my $bit_mis = 0;
	my $found = 0; my $found_peng = 0; 
	
	open my $vcf, '<', $file or die "$!\n";
	while(<$vcf>){
		# make sure its a chromosome
		next if /^#/;
		
		# get required information
		my @fields = split(/\t/,$_);
		my $chr = $fields[0]; $chr =~ s/^chr//;
		my $loc = $fields[1];
		my $key = "$chr:$loc";

		# if we need to guess, extra step
		if($guess_hg){
			my $continue = guessHG('key'=>$key,'ref'=>$fields[3],'bam'=>0);
			
			# determine if change of reference was 
	 		# detected
			if(! $continue){
		 		vcf();
				return;
			}
		}

		# determine if this is marker
		next unless exists $bit_loc{$key};

		my $bit = $bit_loc{$key}-1;
		my $mis = ($fields[$readCol] =~ $misR)? 1:0;
			  
		# determine if sex chromosome marker
		if(exists $allosomes{$key} && $sex){
			delete $ref_allel{$key};

			xmb_builder( 'key' => $key,'data' => $_ );

			next;
		}

		# if its not a sex chromosome marker, then
		# it must be an autosome
		next unless exists $autoSomes{$key};

		my $zyg = ($fields[$readCol] =~ $homo1 or $fields[$readCol] =~ $homo2)? 0:1;

		#fill in the correct amb bit
		$AMB[$bit] = $zyg * (!$mis);
		my $isPengelly = $pengelly{$key};
		$found_peng += $isPengelly;   

		# determine if missed read
		if($mis){
			$missing_markers++;
			$bit_mis = $bit;

			if($isPengelly){$only_peng = 0;}

		}elsif(!$isPengelly){$only_peng = 0;}
		
		$found += (!$mis);
	}

	# if there was a flag to match
	# sex information to the reference genome, then
	# loop through and do so
	if($ref && $sex){
		foreach my $key (keys %ref_allel){
			my $bit = $bit_loc{$key} - 1;
			my ($anc,$alt) = split(":",$allosomes{$key});
			 
			# determine if ref is ancestral or alternate
			# and update the bit accordingly
			 
			# homozygous ancestral
			# check if y is present, if so
			# no modification to the array (already set to 0)
			if($ref_allel{$key} eq $anc && !$SMB[0]){
				$XMB[$bit] = 0;
			}

			# homozygous alternate
			else{
				$XMB[$bit+1] = 1;
			}
		}
	}

	# deal with missing markers
	if(!$ref){
		$missing_markers += 58 - $found;
	}
	$only_peng = ($only_peng && $found_peng == 24);
	
	genMADIB('miss_count'=>$missing_markers,'oPeng'=>$only_peng,'bMis'=>$bit_mis);   
}


=comment

	UTILITY/HELPER SUBROUTINES 
		
		- not directly involved in the process
		- but do help do stuff

	guessHG
		makes an educated guess at the correct human
		genome version, if not specified

		@param key     a string of the form <chr>:<position>
		@param ref     the refrence allel from the sample file entered
		@param bam     1|0 => bam sample | not a bam sample
		@return        1|0 => continue | do not continue

	genMADIB
		generates the MADIB for all forms of input files

		@param miss_count    number of missing bits
		@param oPeng         1|0 only the Pengelly markers read
		@param bMis          the last location of a missed read
		@return              void
	
	xmb_builder
		For vcf and tabix files, this will build the 
		xmb marker bits

		@param key        a string of the form <chr>:<position>
		@param data       VCF row for a specific location
		@return           void

	bam_zygosity
		For bam files, this will automate the quering of 
		the bam file and determine if the specific diploid
		is homozygous or not

		@param chr        the chromosome location
		@param loc        the location on the chromosome
		@param smb        1|0, determine sex or not
		@param anc        ancestral allele
		@param alt        alternate allele
		@return           returns the zygosity of the snp
	
	detReadCol
		For vcf/tbi files, there are more than one samples that can
		be contained. this will determine which coloumn is required 
		to be read

		@param 0				the type of file (tbi or vcf)
		@return 				void
=cut

sub guessHG{
	# declare input variables
	my (%input) = @_;
	my $key = $input{'key'};
	my $smp_ref = $input{'ref'};
	my $isBam = $input{'bam'};
	
	# determine ref allele for initial guess
	my $cur_ref = (exists $ref_allel{$key})? 
		$ref_allel{$key} : "";

	# check reference alleles at hg positions
	# which everone maps out, works
	foreach my $hg (keys %ref_allel_static){
		next if $hg eq $ref_hg;

		my (%ref) = @{$ref_allel_static{$hg}};
		
		# ensure the position exists
		next unless exists $ref{$key};
		my $ref_allel = $ref{$key};

		# determine if ref allel lines up
		if($ref_allel eq $smp_ref){

	 # make sure that initial guess is wrong
	 if($ref_allel eq $cur_ref){
		 return 1;
	 }
	 
	 # if the initial guess is wrong, then
	 # it must be this reference allel
	 $guess_hg = 0;
	 
	 # change to correct reference genome
	 $ref_hg =$hg;
	 (%pengelly) = @{ $pengelly_static{ $hg } };
	 (%bit_loc) = @{ $bit_loc_static{ $hg } };
	 (%autoSomes) = @{ $autoSomes_static{ $hg } };
	 (%ref_allel) = @{ $ref_allel_static{ $hg } };
	 
	 # signal parent subroutine to stop running
	 return 0;
		}
		
	}

	# ensure that the initial guess has an key
	# stored as a location
	if(! exists $bit_loc{$key}){
		return 1;
	}

	# by this point, only the initial guess was
	# correct, thus we are good
	if($cur_ref eq $smp_ref){
		$guess_hg = 0;
		return 1;
	}
	
	# The reference genome could not be determined
	# so die
	die 'Reference genome version not determined, please specify';
}

sub genMADIB{
	my (%input) = @_;
	
	# if no missing markers, MADIB = 0
	if($input{'miss_count'} == 0){
		$MADIB = "000000";
	}
	# if greater than 4, set to 63
	elsif($input{'miss_count'} > 4){
		$MADIB = sprintf("%b",63);
	}
	# only Pengelly, set to 59
	elsif($input{'oPeng'}){
		$MADIB = sprintf("%b",59);
	}
	# only one missing marker
	elsif($input{'miss_count'} == 1){
		$MADIB = "0" x (6-length(sprintf("%b",$input{'bMis'}))) . sprintf("%b",$input{'bMis'});
	}
	# else just show the number of missing markers
	else{
		$MADIB = "0" x (6-length(sprintf("%b",$input{'miss_count'}))) . sprintf("%b",$input{'miss_count'});
	}
}

sub xmb_builder{
	my (%input) = @_;

	# start to define useful constants
	# inorder to build the XMB given
	# a vcf line and allosome key
	my $key = $input{'key'};
	my ($chr,$pos) = split(':',$input{'key'});
	my $row = $input{'data'};
	my @col = split(/\t/,$row);
	my $mis = 0;

	# get information that directly is involved
	# with the value/location of the bit
	my ($anc,$alt) = split(":", $allosomes{$key});

	# determine if it is a Y chromosome
	if($chr eq 'Y'){
		# check to ensure it was not a missed read
		next if $col[$readCol] =~ $misR;
		my($bit,$depth,$marker_id) = split(":", $allosomes{$key});

		# ensure that all trees are the same
		if($Y_marker_name ne 1){
			my ($tree, $branch) = split("-",$marker_id);

			# check to see if tree is the same
			if( index($Y_marker_name, $tree) != -1 || $Y_marker_name eq ""){
				$Y_marker_name = $tree;
			}
			else{
				$Y_marker_name = 1;
			}
		}

		# we only take the deepest branch
		next unless $depth > $max_depth;
		$max_depth = $depth;
		
		$Y_marker_id = $bit;	
		$SMB[0] = 1;

		next;
	}

	my $bit = $bit_loc{$key} -1;

	# determine if multiple alternates were listed
	if(! $col[4] =~ /^\s*\w{1}\s*/){
		$ucn =1;

		# determine if alternate is listed
		if( index($col[4],$alt) != -1){
			$col[4] = $alt;
		}
		else{ $col[4] = "P"; }
	}

	$mis = ($col[$readCol] =~ $misR)? 1:0;

	# homozygous to ancestral allele
	if($col[$readCol] =~ $homo2 && $col[3] eq $anc){
		$XMB[$bit] = 0;
	}

	# homozygous to alternate allele
	elsif($col[$readCol] =~ $homo1 && $col[4] eq $alt){
		$XMB[$bit+1] = 1;
	}

	# heterozygous
	elsif(!$mis && $col[3] eq $anc && $col[4] eq $alt){
		$XMB[$bit+1] = 1;
		$XMB[$bit +1] = 1; $XMB[$bit] = 0;
	}
}

sub bam_zygosity{
	my (%input) = @_;
	
	# get aligned reads at given position
	my @pairs = $file->get_features_by_location(-type=>'read_pair',-seq_id=>$input{'chr'},
						-start=>$input{'loc'},-end=>$input{'loc'});
	# counters
	my $reads = 0;
	my %alleles = ('G' => 0, 'C'=> 0, 'T'=> 0, 'A' => 0);
	
	# loop through pairs
	for my $pair(@pairs){
		my ($f,$s) = $pair->get_SeqFeatures;
		my $start = $input{'loc'} - $f->start;
		
		# filter if score is not good enough
		my $score = ($f->qscore)[$start];
		next unless $score >= $baq;
		
		# increment the allele counters
		my $bp = substr($f->query->dna,$start,1);
		$alleles{$bp}++;
		
		$reads++; #increment reads, accounts for noise
	}
	
	# determine the 2 most populated alleles
	# and compare to see if homo or hetero
	my $al_1 = (sort {$alleles{$a} <=> $alleles{$b}} keys %alleles)[(keys %alleles)-1];
	my $al_2 = (sort {$alleles{$a} <=> $alleles{$b}} keys %alleles)[(keys %alleles)-2];
	my $al1 = $alleles{$al_1};
	my $al2 = $alleles{$al_2};

	# ensure there were reads
	if ( $al1 ==0 && $al2 == 0  ){
		return (exists $input{'smb'})? "1:0" :-1;
	}

	# ensure that it was greater than noise
	my $zyg_prob = $al1/($al1 + $al2) * 100;
	if($zyg_prob**$reads < $noise){
		return (exists $input{'smb'})? "1:0": 0;
	}
		
	# handle if this is sex chromosome bits
	if(exists $input{'smb'} ){
		if($input{'chr'} eq 'Y'){
			if($al_1 eq $input{'anc'} && $al1 > 0.9*($al1 + $al2) ){
				return 1;
			}
			return 0;
		}
	 
		# determine if homozygous
		if($al1 >= 0.90 * ($al1 + $al2)){
			# if homo, determine for which allele
			if($al_1 eq $input{'anc'}){
				return "0:0";
			}
			elsif($al_1 eq $input{'alt'}){
				return "1:1";
			}
			return "1:0";
		}
		
		# then must be some form of 
		# hetero zygous
		else{
			# ensure its hetero between the 
			# ancestral and alternate alleles
			if( $alleles{$input{'anc'}} == 0  or
				$alleles{$input{'alt'}} == 0 ){
					return "1:0";
			}
			
			# heterozygous
			return "0:1";
		}
		 
	}

	# return zygosity
	# ensuring that $al1 is over 90% of reads
	return ( $al1  >= 0.90 * ($al1+$al2) )? 0:1; 
}

sub detReadCol{
	my $type = $_[0];
	my @header;

	if($type eq "tbi"){
		@header = split(/\t/, (split(/\n/,`tabix -H $file`))[-1]);
	}
	else{
		my $header = "";
		open my $vcf, '<', $file or die "$!\n";
		while(<$vcf>){

			# ensure that it is a header
			last unless /^#/;
			$header = $_;
		}
		@header = split(/\t/, $header);
	}
	$readCol = first_index{ /$sampleName/ } @header;
	die "No Sample: $sampleName in specified file\n" unless $readCol > 0;
}


package main;

my $vcfFile = "/export/home/yusuf/geneomeID/sample1.bam";
#   $vcfFile = "/export/home/yusuf/geneomeID/genomeID/samples/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz";
	#$vcfFile = "/export/home/yusuf/geneomeID/HG00157.1000g.vcf.gz";
#	$vcfFile = "/export/home/yusuf/geneomeID/sample1.vcf";

my $genID = genomeID::generate_id('type'=>'bam','file'=>$vcfFile,'hg'=>'hg19','sex'=>1,'ref'=>1);

print "$genID\n";






