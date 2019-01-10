#include <stdio.h>
#include <ctype.h>

const char     *rnames[] =
{
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "UNK", "UNK"
};

/* Convert AA letter to numeric code (0-22) */
int
                aanum(int ch)
{
    const static int aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2,
	22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 22);
}

main()
{
    char ch, id[80], desc[1024], seq[5000];
    int i, nres=0;

    while (!feof(stdin))
    {
	while (getchar() != '\n');
	while ((ch = getchar()) != '>' && ch != EOF)
	    if (isalpha(ch))
	    {
		seq[nres++] = aanum(ch);
	    }
    }

    for (i=0; i<nres; i++)
	printf("%s%s", rnames[seq[i]], i%13 != 12 ? " " : "\n");
}
