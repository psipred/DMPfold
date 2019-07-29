/* pdbhfilt - strip hydrogen atom records from PDB file */

/* David Jones, June 2008 */

/* Copyright (C) 2008 University College London */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Strip hydrogens and waters (& do some other cleanups) from PDB files */

int             main(int argc, char **argv)
{
    int currid=0, atomid=1;
    char buf[256], lastid[10], idstr[10];

    if (argc > 1)
	currid = atoi(argv[1])-1;

    while (!feof(stdin))
    {
	if (!fgets(buf, 256, stdin))
	    break;
	if (!strncmp(buf, "ATOM ", 5) || !strncmp(buf, "HETATM ", 7))
	{
	    if (!strncmp(buf+17, "HOH", 3) || !strncmp(buf+17, "WAT", 3) || !strncmp(buf+17, "DOD", 3))
		continue;
	    strcpy(buf+54, "  1.00  0.00\n");
	    if (buf[12] == 'H' || buf[13] == 'H')
		continue;
	    if (buf[12] == 'D' || buf[13] == 'D')
		continue;
	    if (!strncmp(buf+13, "CD  ILE", 7))
		strncpy(buf+13, "CD1 ILE", 7);
	    if (!strncmp(buf+13, "O1 ", 3))
		strncpy(buf+13, "O  ", 3);
	    if (!strncmp(buf+13, "O2 ", 3))
		continue;
	    if (strncmp(buf+21, lastid, 5))
	    {
		strncpy(lastid, buf+21, 5);
		currid++;
	    }
	    if (argc > 1)
	    {
		sprintf(idstr, "%5d", currid);
		strncpy(buf+21, idstr, 5);
	    }
	    sprintf(idstr, "%5d", atomid++);
	    strncpy(buf+6, idstr, 5);
	}
	printf("%s", buf);
    }
    return 0;
}
	
