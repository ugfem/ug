/****************************************************************************/
/*                                                                          */
/* File:      MacGui.m                                                      */
/*                                                                          */
/* Purpose:   definition of constants used in resource file                	*/
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            6900 Heidelberg                                               */
/*            internet: bastian@iwr1.iwr.uni-heidelberg.de                  */
/*                                                                          */
/* History:   27.02.92 begin, ug version 2.0                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/* ids for 'WIND' resources */
#define appWinId		129
#define aboutWinId		130

/* ids for 'CNTL' resources */
#define vScrollBarId	128
#define hScrollBarId	129

/* ids for pltt resources */
#define defaultPaletteId	1000

/* ids for PICT resources */
#define kawPict				128
#define addressPict			129
#define ncsaPict			130
#define toolboxPict			150

/* ids for CURS resources */
#define watchCurs			128
#define textCurs			129
#define	arrowCurs			150
#define crossCurs			151
#define blobCurs			152 
#define handCurs			153
#define triangleCurs		154

/* ids for DLOGs */
#define newDialog			128
#define openDialog			129
#define threeButtonDialog   130
#define twoButtonDialog     131
#define oneButtonDialog     132
#define saveOptionsDialog	133
#define viewOptionsDialog	134
#define plotOptionsDialog	135
#define factorDialog		136
#define positionDialog		137
#define innerDialog			138
#define selectSegmentDialog 139
#define moveBndNodeDialog	140
#define moveInnerNodeDialog 141
#define uiSettingsDialog	142

/* ids for DITLs */

/* general : */
#define OKButton 			1
#define cancelButton		2
#define YesButton			1
#define NoButton			3

/* DITL id 128/129 (newDialog/openDialog) */
#define newdocTextItem		3
#define heapSizeItem		7
#define problemItem			8
#define formatItem			9
#define dataButtonItem		10

/* DITL id 133 */
#define saveDataButton		4
#define saveFlatButton		5
#define commentText			7

/* DITL id 134 */
#define nodesItem			4
#define edgesItem			5
#define elementsItem		6
#define boundaryItem		7
#define idsItem				8
#define marksItem			9
#define useColorItem		10

/* DITL id 135 */
#define saveContoursButton	3
#define loadContoursButton	4
#define plotProcPopUp		5
#define colorFromText		6
#define colorToText			7
#define depthText			8
#define contoursOnRadio		9
#define numOfLinesText		10
#define equidistantRadio	11
#define equiFromText		12
#define equiToText			15
#define customRadio			13
#define lineNoText			14
#define valueText			16
#define imageOnRadio		17
#define imageFromText		18
#define imageToText			19
#define RangeButton			33

/* DITL id 136 */
#define factorHeading		3
#define factorText			4

/* DITL id 137 */
#define positionPopUp		4

/* DITL id 138 */
#define innerPosItem		4

/* DITL id 139 */
#define selSegText			2
#define selSegPopUp			3

/* DITL id 140 */
#define segPopUp			4
#define parameterText		6

/* DITL id 141 */
#define xTextItem			5
#define yTextItem			7

/* DITL id 142 */
#define clicktolTextItem	5
#define bndresTextItem		7

/* MENU resources */
#define menuCount			3

/* menu ids */
#define	appleID		128
#define	fileID		129
#define	displayID	130

#define	appleM		0
#define fileM		1
#define	displayM	2

/* apple menu commands */
#define aboutCommand		1

/* file menu commands */
#define openCommand			1
#define closeCommand		2
#define quitCommand			4

/* display menu commands */
#define refreshCommand		1
#define hdfCommand			2
#define pictCommand			3
