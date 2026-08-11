/************************************************************************/
/*  iso532_oracle - reference-data generator for the SQAT test suite.    */
/*                                                                      */
/*  Thin driver around the UNMODIFIED ISO 532-1:2017 Annex A.4          */
/*  reference library (ISO_532-1.c / ISO_532-1.h). It calibrates an     */
/*  input WAV exactly as SQAT's calibrate.m does, calls the reference   */
/*  implementation, and writes the result as CSV for use as golden      */
/*  data by the MATLAB tests.                                          */
/*                                                                      */
/*  The ISO reference sources are NOT vendored here: they ship with the */
/*  standard and are not ours to redistribute. Point ISO_SRC at your    */
/*  copy to rebuild the golden data (see Makefile). The generated CSVs  */
/*  in test/golden/ are committed so the test suite runs without them.  */
/*                                                                      */
/*  Usage:                                                              */
/*    iso532_oracle tv <sig.wav> <cal.wav> <cal_dB> <out.csv> [spec.csv]*/
/*        time-varying loudness, N(t) at SR_LEVEL = 2000 Hz            */
/*    iso532_oracle st <sig.wav> <cal.wav> <cal_dB> <skip> <out.csv>    */
/*        stationary loudness from a signal, single value              */
/*    iso532_oracle lv <levels.txt> <out.csv> [spec.csv]                */
/*        stationary loudness from 28 third-octave levels              */
/*                                                                      */
/*  Add D as a trailing argument to any subcommand for diffuse field.   */
/************************************************************************/

#include "ISO_532-1.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*  Rate to which the reference main program downsamples the summed loudness
    for output. Defined in the ISO Annex A.4 ISO_532-1_main.c, not in the
    library header, so it is repeated here.                                */
#define SR_LOUDNESS 500

/*  ---------------------------------------------------------------- */
/*  Minimal WAV reader: mono/stereo, 16-bit PCM or 32-bit IEEE float. */
/*  Stereo takes channel 1, matching the standard's "one channel"    */
/*  requirement and audioread(...,1) on the MATLAB side.             */
/*  ---------------------------------------------------------------- */
static int rd_wav(const char *path, double **out, int *n, double *sr)
{
    FILE *f = fopen(path, "rb");
    unsigned char h[8], b[16];
    char id[5] = {0};
    int bits = 0, ch = 0, i;
    long rate = 0;

    if (!f) { fprintf(stderr, "iso532_oracle: cannot open %s\n", path); return -1; }
    if (fseek(f, 12, SEEK_SET) != 0) { fclose(f); return -1; }   /* RIFF....WAVE */

    for (;;)
    {
        unsigned int sz;
        if (fread(h, 1, 8, f) != 8) { fclose(f); return -1; }
        memcpy(id, h, 4);
        sz = h[4] | (h[5] << 8) | (h[6] << 16) | ((unsigned)h[7] << 24);

        if (!strcmp(id, "fmt "))
        {
            if (fread(b, 1, 16, f) != 16) { fclose(f); return -1; }
            ch   =  b[2] | (b[3] << 8);
            rate =  b[4] | (b[5] << 8) | (b[6] << 16) | ((long)b[7] << 24);
            bits = b[14] | (b[15] << 8);
            if (sz > 16) fseek(f, (long)sz - 16, SEEK_CUR);
        }
        else if (!strcmp(id, "data"))
        {
            int ns = (int)(sz / (unsigned)(bits / 8) / (unsigned)ch);
            double *d = (double*)malloc(sizeof(double) * (size_t)ns);
            if (!d) { fclose(f); return -1; }
            for (i = 0; i < ns; i++)
            {
                if (bits == 16) { short v; if (fread(&v,2,1,f)!=1) break; d[i] = v / 32768.0; }
                else            { float v; if (fread(&v,4,1,f)!=1) break; d[i] = v; }
                if (ch == 2) fseek(f, bits / 8, SEEK_CUR);   /* skip channel 2 */
            }
            fclose(f);
            *out = d; *n = ns; *sr = (double)rate;
            return 0;
        }
        else fseek(f, (long)sz + (sz & 1), SEEK_CUR);
    }
}

/*  Calibration factor, identical to SQAT utilities/calibrate.m:
    CalFactor = sqrt( 10^(L/10) * I_REF / mean(ref^2) )                */
static double cal_factor(const double *ref, int n, double level_dB)
{
    double ms = 0; int i;
    for (i = 0; i < n; i++) ms += ref[i] * ref[i];
    return sqrt(pow(10.0, level_dB / 10.0) * I_REF / (ms / n));
}

/*  Read 28 third-octave levels; "<f> : <level>" per line, # = comment  */
static int rd_levels(const char *path, double *lvl)
{
    FILE *f = fopen(path, "r");
    char line[256];
    int k = 0;
    if (!f) { fprintf(stderr, "iso532_oracle: cannot open %s\n", path); return -1; }
    while (fgets(line, sizeof line, f) && k < N_LEVEL_BANDS)
    {
        char *c = strchr(line, ':');
        if (line[0] == '#' || !c) continue;
        lvl[k++] = atof(c + 1);
    }
    fclose(f);
    if (k != N_LEVEL_BANDS)
    {
        fprintf(stderr, "iso532_oracle: expected %d levels, got %d\n", N_LEVEL_BANDS, k);
        return -1;
    }
    return 0;
}

static void wr_spec(const char *path, double *S[N_BARK_BANDS], int nt)
{
    FILE *o; int b, t;
    if (!path) return;
    o = fopen(path, "w");
    fprintf(o, "bark");
    for (t = 0; t < nt; t++) fprintf(o, ",t%d", t);
    fprintf(o, "\n");
    for (b = 0; b < N_BARK_BANDS; b++)
    {
        fprintf(o, "%.1f", (b + 1) / 10.0);
        for (t = 0; t < nt; t++) fprintf(o, ",%.10f", S[b][t]);
        fprintf(o, "\n");
    }
    fclose(o);
}

int main(int argc, char **argv)
{
    double *N = NULL, *S[N_BARK_BANDS];
    int i, nlev = 1, ret, field = SoundFieldFree;
    const char *mode, *specout = NULL, *out = NULL;

    if (argc < 3) { fprintf(stderr, "usage: see header of %s\n", __FILE__); return 2; }
    mode = argv[1];
    for (i = 1; i < argc; i++) if (!strcmp(argv[i], "D")) field = SoundFieldDiffuse;

    if (!strcmp(mode, "lv"))                        /* ---- from levels ---- */
    {
        double lvl[N_LEVEL_BANDS], *TOL[N_LEVEL_BANDS];
        if (argc < 4) { fprintf(stderr, "usage: lv <levels.txt> <out.csv> [spec.csv]\n"); return 2; }
        if (rd_levels(argv[2], lvl)) return 1;
        out = argv[3]; specout = (argc > 4 && strcmp(argv[4],"D")) ? argv[4] : NULL;

        for (i = 0; i < N_LEVEL_BANDS; i++) TOL[i] = &lvl[i];
        N = (double*)calloc(1, sizeof(double));
        for (i = 0; i < N_BARK_BANDS; i++) S[i] = (double*)calloc(1, sizeof(double));
        ret = f_loudness_from_levels(TOL, 1, field, LoudnessMethodStationary, N, S);
        nlev = 1;
    }
    else                                            /* ---- from signal ---- */
    {
        double *sig, *cal, sr, csr, k, skip = 0;
        int ns, nc, meth;
        struct InputData in;

        if (argc < 6) { fprintf(stderr, "usage: see header of %s\n", __FILE__); return 2; }
        if (rd_wav(argv[2], &sig, &ns, &sr))  return 1;
        if (rd_wav(argv[3], &cal, &nc, &csr)) return 1;

        k = cal_factor(cal, nc, atof(argv[4]));
        for (i = 0; i < ns; i++) sig[i] *= k;

        in.NumSamples = ns; in.SampleRate = sr; in.pData = sig;

        if (!strcmp(mode, "st"))
        {
            meth = LoudnessMethodStationary;
            skip = atof(argv[5]);
            out  = argv[6];
            nlev = 1;
        }
        else
        {
            meth = LoudnessMethodTimeVarying;
            out  = argv[5];
            specout = (argc > 6 && strcmp(argv[6],"D")) ? argv[6] : NULL;
            nlev = ns / (int)(sr / SR_LEVEL);
        }

        N = (double*)calloc((size_t)nlev, sizeof(double));
        for (i = 0; i < N_BARK_BANDS; i++) S[i] = (double*)calloc((size_t)nlev, sizeof(double));

        ret = f_loudness_from_signal(&in, field, meth, skip, N, S, nlev);
    }

    if (ret < 0) { fprintf(stderr, "iso532_oracle: reference library error %d\n", ret); return 1; }

    /*  Time-varying results are emitted on the SR_LOUDNESS = 500 Hz grid,
        which is what the ISO reference main program reports and what SQAT's
        OUT.time / OUT.InstantaneousLoudness use. Stationary results are a
        single value and are written as-is.                                 */
    {
        int dec = (ret > 1) ? (SR_LEVEL / SR_LOUDNESS) : 1;
        int nw  = 0;
        FILE *o = fopen(out, "w");
        if (!o) { fprintf(stderr, "iso532_oracle: cannot write %s\n", out); return 1; }
        fprintf(o, "t,N\n");
        for (i = 0; i < ret; i += dec, nw++)
            fprintf(o, "%.4f,%.8f\n", i / (double)SR_LEVEL, N[i]);
        fclose(o);
        fprintf(stderr, "iso532_oracle: %s -> %s (%d sample%s @ %d Hz)\n",
                mode, out, nw, nw == 1 ? "" : "s",
                (ret > 1) ? SR_LOUDNESS : 1);
    }
    wr_spec(specout, S, ret);
    return 0;
}
