/*

	Simulatore di scontri galattici.
	Evoluzione di un antico (e dico antico) algoritmo in GW-Basic,
	divenuto famoso ma scritto non-so-da-chi.
	A parte l'idea, comunque, StarWays non ha nulla di quell'algoritmo,
	neanche il nucleo di calcolo, che qui ä ottimizzato e tradotto in
	assembly per la fpu.

	Il concetto ä lo stesso, ma StarWays ä un'altra cosa.
	Su "la vetrina del software" di Nuovo Orione, recensirono entusiasti
	il programma in basic suo avo ancestrale, ed io ho scritto StarWays
	per farmi pubblicitÖ ed ammaliare i recensori di Nuovo Orione:
	se si sono divertiti con il Basic, figuriamoci come accetteranno
	StarWays; oltre a questo, perï, c'ä stato il mio interesse personale
	nel voler avere, per primo, un "giocattolo tecnologico" come questo.

	Studiato come pubblico dominio, anche perchÇ non mi ha impegnato
	poi molto, calcoliamo una settimanetta. Perï forse conviene
	registrarlo, perchÇ... hai visto mai? Se piacesse molto, lo potrei
	triplicare a livello di funzioni e di sofisticatezza, per poi
	rilasciare una versione commerciale. Oh, yeah. Non male.

	Tutto fatto.

	13 febbraio 2003: tradotto in inglese (solo stringhe).
	february 13th, 2003: translated to english (strings only).

*/

#include <assembly.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <alloc.h>
#include <math.h>
#include <time.h>
#include <dir.h>
#include <dos.h>
#include <io.h>

double mindiff = 1;
unsigned long global_x, global_y;
long stk_lbx = -319, stk_lby = -239, stk_ubx = 319, stk_uby = 239;

/* Header per files di tipo BMP. */

char bmpheader[54] = { 0x42, 0x4d, 0x76, 0x58, 0x02, 0x00, 0x00, 0x00,
		       0x00, 0x00, 0x76, 0x00, 0x00, 0x00, 0x28, 0x00,
		       0x00, 0x00, 0x80, 0x02, 0x00, 0x00, 0xE0, 0x01,
		       0x00, 0x00, 0x01, 0x00, 0x04, 0x00, 0x00, 0x00,
		       0x00, 0x00, 0x00, 0x58, 0x02, 0x00, 0x13, 0x0b,
		       0x00, 0x00, 0x13, 0x0b, 0x00, 0x00, 0x00, 0x00,
		       0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

/* Limiti. */

#define M_GAL	    4 // Massimo numero di galassie ammesse.

/* Definizioni utili. */

#define Ex	    0
#define Sa	    1
#define Sc	    2
#define S4	    3

char *gtyp[] = { "EX ELLIPTICAL    ",
		 "SA 2-ARMS SPIRAL ",
		 "SC THINNER SPIRAL",
		 "S4 4-ARMS SPIRAL " };

#define NORMALE     0 // Per la funzione "Testo".
#define CENTRATO    3
#define ALL_DX      4

#define DEG  180/M_PI // Per convertire radianti in gradi.

/* Variabili fondamentali e costanti fisiche necessarie.
   Comprensive dei valori iniziali. */

unsigned seed = 4112;	         // Seme random per il generatore pseudo.

unsigned g  = 2;    	         // Numero di galassie.
unsigned ls = 6000;	         // Limite nr. stelle = 48000 / 8 = 6000.
unsigned ms = ls / M_GAL;        // Massimo stelle per galassia.

const double MIL = 3.1536e+13;   // Secondi in un milione di anni.

double single_step = MIL;        /* Secondi tra un fotogramma e l'altro.
				    (Default un milione di anni). */

const double C = 2.99574e+5;	 // VelocitÖ della luce in km/sec.

const double G = 6.67e-11;	 // Costante di gravitazione universale.

const double A = 9.4608e+15;	 // 1 anno-luce espresso in metri.

const double M = 1.9701e+30;     // 1 massa solare espressa in kg.

double S = 1000 * single_step;   /* Migliaia di secondi reali in un passo
				    d'animazione. Dividendo le velocitÖ
				    effettive (espresse in Metri per Step,
				    che divisi per 1000 danno i Km per Step,
				    che divisi ancora per single_step danno
				    il risultato finale) per questo valore,
				    si ottiene il corrispondente valore
				    espresso in Km/Sec, ovvero in una scala
				    pió facilmente comprensibile.
				    PerchÇ "al secondo"? In effetti, il
				    valore realmente visibile dipende da
				    quanti steps di calcolo vengono svolti
				    in un secondo reale: questo numero (pcs)
				    ä visualizzato dal programma di continuo,
				    percui ci se ne puï fare un'idea.
				    D'altronde, per poter studiare il moto
				    stellare si deve leggere la velocitÖ che
				    avrebbero le stelle in una situazione
				    reale: l'accelerazione del tempo nella
				    simulazione serve solo per comoditÖ.
				    E` poco comprensibile e poco pratico
				    far avanzare la simulazione esattamente
				    ad un fotogramma (o step) al secondo,
				    mentre suddividere single_step in
				    fammenti dati da quanti fotogrammi
				    vengono visualizzati in un secondo non
				    ä affidabile (perchÇ il timer interno
				    ä piuttosto impreciso e perchÇ in effetti
				    il valore di pcs dipende dalle condizioni
				    nelle quali si assiste alla simulazione.
				    Pió stelle sullo schermo rallenteranno
				    tutto, al pari dell'attivazione di una
				    delle tre funzioni di inseguimento. */

double TT = 0; /* Tempo trascorso nella simulazione (in secondi). */

/* Dati sulle galassie. */

unsigned n[M_GAL] = {   ms,  ms/5 }; // Numero di stelle per galassia.
unsigned t[M_GAL] = {   Sa,    Ex }; // Archetipo della galassia.
unsigned X[M_GAL] = { 4112,  4112 }; // Seme della galassia.
double   d[M_GAL] = {  5E4,   5E4 }; // Diametro galassia.
double   s[M_GAL] = {  5E8,   1E8 }; /* Massa dell'attrattore centrale.
				      Nota: Per default, si considera un buco
					    nero ipermassivo di massa pari a
					    mezzo miliardo di Soli. */

int      trgt[M_GAL] = { 0, 0 }; // Eventuale bersaglio.
double   ashf[M_GAL] = { 0, 0 }; // Alfa-Shift (correzione su declinazione).
double   bshf[M_GAL] = { 0, 0 }; // Beta-Shift (correzione su ascensione retta).
double   vlct[M_GAL] = { 0, 50*S }; // VelocitÖ del nucleo verso il bersaglio. (2500*S x save me).

double   incl[M_GAL] = { 0, 0 }; // inclinazione della galassia.
double   tilt[M_GAL] = { 0, 0 }; // rotazione della galassia sul suo asse.
char  	 acws[M_GAL] = { 0, 0 }; // flag: rotazione stelle antioraria.

double gx[M_GAL]  = { 0 * A, 2E4 * A }, // 2nd: 0 x save me.
       gy[M_GAL]  = { 0 * A, 1E5 * A }, // 2nd: -3.15E6 x save me.
       gz[M_GAL]  = { 0 * A, 0 * A }; // Posizione (iniziale).
double gvx[M_GAL] = { 0 * A },
       gvy[M_GAL] = { 0 * A },
       gvz[M_GAL] = { 0 * A }; // VelocitÖ proiettata sugli assi X, Y, Z.
double bgx[M_GAL] = { 0 * A, 2E4 * A },
       bgy[M_GAL] = { 0 * A, 1E5 * A },
       bgz[M_GAL] = { 0 * A, 0 * A }; // Copia di backup delle posizioni.

/* Dati sui punti massa, ovvero sulle stelle.
   La massa (m) decide tutto: raggio, classe ed emissione luminosa. */

float  far *m;	 	  // Massa (0.1 - 1000) espressa  in masse solari.
double far *x,  *y,  *z;  // Posizione (coordinate tre-d) in anni-luce.
double far *vx, *vy, *vz; // VelocitÖ (componenti tre-d)  in anni-luce / step.

// Variabili di ottimizzazione (per i cicli assembly di calcolo).

unsigned seg_adapted, seg_m;
unsigned seg_x, seg_y, seg_z;
unsigned seg_vx, seg_vy, seg_vz;

/* Aree di supporto funzioni mem & rcl. */

double mcam_x, mcam_y, mcam_z;
double m_alfa, m_beta;

/* Dati sul punto d'osservazione e sulla proiezione 3d-2d. */

double cam_x = 0;           // Posizione videocamera virtuale.
double cam_y = -40000 * A;
double cam_z = 0;

double alfa  = M_PI_2;      // Angolazione.
double beta  = 0;

double delta_alfa = 0;      // Variazione angolazione.
double delta_beta = 0;

double _alfa, _beta;	    // Angolazione salvata in un backup.

double delta_pcam = 0;      // Accelerazione dell'osservatore.

double ob = 300; /* Obiettivo usato, o distanza dal piano di proiezione.
		    - A valori pió piccoli corrispondono "grandangolari".
		    - A valori pió grandi corrispondono "teleobiettivi".
		    - Al valore "300" corrisponde un normale 50 mm. */

double min_dist = ob * A; // Distanza normale minima di una stella visibile.

double cosbeta, sinbeta; // Variabili temporanee di ottimizzazione.
double cosalfa, sinalfa; // XchÇ le funzioni trigonometriche sono tartarughe.

double cosbeta_x_ob, sinbeta_x_ob; // Come sopra giÖ moltiplicate per la
double cosalfa_x_ob, sinalfa_x_ob; // distanza dal piano di proiezione. OlÇ.

unsigned anchor = 1639; /* Indice della stella alla quale si ä ancorati.
			   -1, ovvero 65535, significa NON ANCORATI. */

char follow = 2; /* ModalitÖ di inseguimento del moto: 0 = nessuna,
		    +1 = inquadratura a prua, -1 = inquadratura a poppa,
		    +2 = inseguimento del nucleo galattico. */

double see_x; /* Coordinate del punto da inquadrare. Funzione "watch". */
double see_y;
double see_z;

double cam_vel    = 0; // VelocitÖ attuale dell'osservatore.
double delta_xpos = 0; // Componente x della velocitÖ.
double delta_ypos = 0; // Componente y della velocitÖ.
double delta_zpos = 0; // Componente z della velocitÖ.
double rif_x	  = 0; // Componente x di riferimento.
double rif_y	  = 0; // Componente y di riferimento.
double rif_z	  = 0; // Componente z di riferimento.

char   wait       = 0; // Segnale per la sospensione della simulazione.

char   help       = 0; // Segnale della presenza delle istruzioni on-line.

char   started	  = 1; // Segnale: la simulazione iniziale non ä modificata.

char   starting   = 1; // Segnale: la simulazione iniziale ä iniziata.

char   grid       = 0; // Segnale: visualizzare la griglia 3d?

char   glow	  = 0; // Segnale: visualizzazione in VESA?

char   hylit	  = 0; // Segnale: VESA in alta luminositÖ con i "trails"?

/* Variabili globali specifiche. */

unsigned char key;  // Tasto premuto nel ciclo principale.
int	      pcs;  // Passi di calcolo per secondo.
clock_t	      tci;  // Variabile per il cronometraggio dei pcs.
unsigned      P=0;  // Precedente condizione dei pulsanti del mouse.

char	      pcsline[50] = "TIMING IN PROGRESS...\0";

double	      fpcs = 30;

double	      bbangle = 0; // Angolo per la rotazione dei bordi delle bolle.

char	      exitflag = 0; // Messaggio d'uscita (richiesta conferma).
char	      msg = 0; // Flag della bolla "Modifica Schemi Galattici".
unsigned      sglx = 0; // Galassia attualmente modificabile.
char          glxpar = 0; // Parametro galattico da modificare.
char	      subglxpar = 0; // Sottoparametro galattico.
char	      ccd_settings = 0; // Impostazioni preparatorie immagini CCD.
char	      ccdpar = 0; // Parametro del CCD da modificare.
char	      lens = 0; // Impostazioni della videocamera: obiettivo.
char	      gas = 0; /* Non Grandi Avventure Sanitarie (torna), e neanche
			  Grandi Avventure A Fumetti (non torna), bensç
			  Gestione Archivio Simulazioni: flag-bolla. */

char	      input[18], ic = 0; // Stringa valori d'ingresso e suo contatore.
char	      sign, point; // Il parametro supporta segni e/o punti decimali?
char	      ics = 0; // Numero di cifre supportate dalla stringa d'ingresso.

int	      errorflag = -1; // Galassia che eventualmente ha superato C.

/* Variabili globali ausiliarie, ad uso generico. */

unsigned      aux;  // Contatore ausiliario interno.
unsigned      pan;  // Contatore ausiliario esterno.

unsigned      i, j; // Indici generici.

unsigned char gby;  // Byte generico.

double        k, q, w, ww; // Variabili generiche a doppia precisione/80 bit.
double	      f1, f2, f3; // Variabili generiche a doppia precisione/64 bit.
double	      f4, f5; // Ausiliarie.

char	      ccd_outfile[13] = "FILENAME";
char	      simfile[13] = "FILENAME";

char	      escfunc[76]; // Buffer testuale: descriz. funzione tasto ESC.
char 	      gpal[768]; // Tavola colori per il modo vesa "glowstar.tab"
char	      tx[107]; // Buffer testuale general-purpose.

/* Dati per tracciare le selle semplici. */

double _m, _x, _y, _z, one = 1; // Variabili di passaggio.
				// "one" ä praticamente una costante,
				// ma come variabile funziona pió velocemente.

double size_a = 0.1;
double size_b = 0.3;
double size_c = 0.7;

/* Dati sulle glowin' stars. */

unsigned gws_szk[13] = {  6, 10, 13, 15, 17, 19, 21,
			 23, 25, 27, 29, 32, 35  };

unsigned gws_offset[13] = { 6558, 6458, 6289, 6064, 5775, 5414,
			    4973, 4444, 3819, 3090, 2249, 1225, 0 };

int spx, spy;          // posizione del centro.
unsigned cr_sp_page;   // pagina in "desperation".

char gws[6595];        // Immagini-base delle glowin' stars: "glowstar.pic"

double size_g0 = 0.02; // Limite inferiore profonditÖ.
double size_g1 = 12;   // Limite superiore profonditÖ.

/* Dati per le stelle CCD. */

double mag_factor = 0;
unsigned exposure = 8;

int pcrlim = 0;
double crlim;

unsigned char far *desperation; // Display 640x480 tutto in convenzionale (!)
unsigned seg_desperation; // suo segmento.

/* Bitmapped 6x6 font (106 col. per 56 righe), monocromatico, per i testi. */

unsigned char font[65*6] = {
 0,  0,  0,  0,  0,  0, // Spazio.
 8,  8,  8,  8,  0,  8, // Punto esclamativo.
20, 20,  0,  0,  0,  0, // Virgolette.
 0, 20, 62, 20, 62, 20, // Cancelletto.
28, 40, 28, 10, 28,  8, // Dollaro.
36,  4,  8,  8, 16, 18, // Percento.
62, 62, 62, 62, 62, 62, // E commerciale (non c'ä: al suo posto c'ä un blocco di colore che uso come cursore).
 4,  8,  0,  0,  0,  0, // Apice.
 4,  8,  8,  8,  8,  4, // Parentesi tonda aperta.
16,  8,  8,  8,  8, 16, // Parentesi tonda chiusa.
40, 16, 40,  0,  0,  0, // Asterisco.
 8,  8, 62,  8,  8,  0, // Pió.
 0,  0,  0,  0,  4,  8, // Virgola.
 0,  0, 62,  0,  0,  0, // Meno.
 0,  0,  0,  0,  0,  8, // Punto.
 4,  4,  8,  8, 16, 16, // Barra.
62, 34, 34, 34, 34, 62, // Zero.
 4, 12,  4,  4,  4, 14, // Uno.
28, 34,  4,  8, 16, 62, // Due.
28,  2, 28,  2, 34, 28, // Tre.
 2,  6, 10, 18, 63,  2, // Quattro.
30, 16, 28,  2,  2, 60, // Cinque.
14, 16, 60, 34, 34, 28, // Sei.
62, 34,  4,  4,  8,  8, // Sette.
28, 34, 28, 34, 34, 28, // Otto.
28, 34, 34, 30,  4, 56, // Nove.
 0,  8,  0,  0,  8,  0, // Duepunti.
 0,  8,  0,  8,  8, 16, // Punto e virgola.
 0,  4,  8, 16,  8,  4, // Minore.
 0,  0, 30,  0, 30,  0, // Uguale.
 0, 16,  8,  4,  8, 16, // Maggiore.
28, 34,  2,  4,  0,  4, // Punto interrogativo.
 0,  0,  8,  8,  0,  0, // Chiocciola (puntino di separazione).
28, 34, 34, 62, 34, 34, // A (uppercase).
60, 34, 60, 34, 34, 60, // B.
28, 34, 32, 32, 34, 28, // C.
60, 34, 34, 34, 34, 60, // D.
62, 32, 48, 32, 32, 62, // E.
62, 32, 48, 32, 32, 32, // F.
28, 34, 32, 38, 34, 28, // G.
34, 34, 62, 34, 34, 34, // H.
 8,  8,  8,  8,  8,  8, // I.
 6,  2,  2,  2, 34, 28, // J.
34, 36, 56, 40, 36, 34, // K.
32, 32, 32, 32, 32, 62, // L.
34, 54, 42, 34, 34, 34, // M.
34, 50, 42, 38, 34, 34, // N.
28, 34, 34, 34, 34, 28, // O.
60, 34, 34, 60, 32, 32, // P.
28, 34, 34, 34, 36, 26, // Q.
60, 34, 34, 60, 34, 34, // R.
28, 32, 28,  2,  2, 28, // S.
62,  8,  8,  8,  8,  8, // T.
34, 34, 34, 34, 34, 28, // U.
34, 34, 34, 20, 20,  8, // V.
34, 34, 42, 42, 54, 34, // W.
34, 20,  8, 20, 34, 34, // X.
34, 34, 20,  8,  8,  8, // Y.
62,  2,  4,  8, 16, 62, // Z.
 0,  0,  0,  0,  0,  0, // Parentesi quadra aperta (manca).
 7,  5,  7,  0,  0,  0, // Barra rovescia (ä lo zerino).
 0,  0,  0,  0,  0,  0, // Parentesi quadra chiusa (manca).
14,  2,  4,  7,  0,  0, // Ordinale femminile (ä l'esponente 2).
 0,  0,  0,  0,  0, 63, // Sottolinea.
 8,  4,  0,  0,  0,  0  // Accento.
};

// Disattivazione dei messaggi DOS sugli errori critici.

#pragma warn -par
int handler (int errval, int ax, int bp, int si)
{
	hardresume (_HARDERR_IGNORE);
	return (_HARDERR_IGNORE);
}
#pragma warn +par

/* Funzioni per l'introduzione, eventualmente in vesa. */

#define _640	0x69
#define _1024	0x6d

unsigned p_align;
unsigned l_align;

unsigned long screenlength;
unsigned long screenheight;

void initvesa (char pixels_on_line)
	/* Apre la modalitÖ vesa indicata con "pixels_on_line":
	   un valore che, preceduto da un sottolinea per definizione,
	   indica da quanti pixels ä composta ogni riga sullo schermo. */
	/* I valori validi sono: _1024 che corrisponde a 1024x768;
				 _640  che corrisponde a 640x480. */
{
	asm {   mov al, pixels_on_line
		xor ah, ah
		int 0x10 }

	switch (pixels_on_line) {
		case _1024:	screenlength = 1024;
				screenheight = 768;
				p_align      = 64;
				l_align	     = 0;
				break;
		case _640:	screenlength = 640;
				screenheight = 480;
				p_align	     = 102;
				l_align	     = 256;
				break;
	}
}

void pselect (char videopage) // Seleziona la pagina indicata.
{
	asm {   mov dl, videopage
		mov ax, 0x4f05
		xor dh, dh
		mov bx, 1
		int 0x10 }
}

void flr_on_vesa (char *nomefile_flr, int x, int y, int l, int h, long offs)
{
	int fh = _open (nomefile_flr, 0);
	int scanline = y;

	unsigned vptr, r, pg;

	if (fh==-1) return;

	if (offs) lseek (fh, offs, SEEK_SET);

	pg   = scanline / p_align;
	r    = pg;

	vptr = (scanline % p_align) * screenlength + x;

	if (vptr < r * l_align) pg--;
	vptr -= r * l_align;

	pselect (pg);

	while (scanline<y+h) {
		_read (fh, tmppal, l);
		asm {   lea si, tmppal
			les di, adaptor
			add di, vptr
			mov cx, l }
	l_draw:	asm {	mov al, [si]
			mov es:[di], al
			inc si
			add di, 1
			jnc cross
			inc pg
			mov ax, 0x4f05
			mov dx, pg
			mov bx, 1
			int 0x10 }
	cross:	asm {	dec cx
			jnz l_draw
			mov vptr, di
			mov ax, word ptr screenlength
			sub ax, l
			add vptr, ax
			jnc jmpr
			inc pg
			mov ax, 0x4f05
			mov dx, pg
			mov bx, 1
			int 0x10 }
	jmpr:	scanline++;
	}

	_close (fh);
}

// Attende un tasto.

void simwait ()
{
	pcs = fpcs;

	while (!tasto_premuto());
	while (tasto_premuto()) attendi_pressione_tasto();
}

/* Funzioni di utilitÖ grafica. */

double   rx, ry, rz, z2, z2d; // Variabili di conversione 3d-2d.
int      x2d, y2d;	      // Coordinate convertite in 2d.
unsigned shift;     	      // Variabile per la funzione "plot".
char	 airbrush = 1;	      // Pennello ad aria per i punti delle linee.

/* Tracciamento bastoncini (parte 2d). */

void Stick (unsigned long x, unsigned long y,
	    unsigned long x2, unsigned long y2)
{
	long a, b, L, temp;

	if (x==x2) x2++;

	a = x2-x;
	if (a<0) {
		temp = x2; x2 = x; x = temp;
		temp = y2; y2 = y; y = temp;
		a = -a;
	}

	L = a;
	b = y2-y;
	if (b>0) {
		if (b>L)
			L=b;
	}
	else {
		if (-b>L)
			L=-b;
	}
	L++;

	x2 <<= 16;
	global_x = x << 16;

	a <<= 16; a /= L;
	b *= 65536; b /= L;
	global_y = y << 16;

	a *= airbrush;
	b *= airbrush;

	/* Traccia il segmento tramite DMA. */

	asm lea si, global_x[2];
	asm lea di, global_y[2];
	asm les ax, dword ptr adapted;

	/*

		Nota importante:

			Tutte le operazioni precedute dal prefisso 0x66
			vanno considerate con operandi di dimensioni doppie.
			Difatti, poichÇ l'assemblatore in-line non supporta
			i registri a 32 bit, ä necessario aggiungere i
			prefissi d'estensione manualmente.
			Cosç facendo, i puntatori di tipo word (word ptr)
			vanno considerati come dword ptr, mentre i registri
			sono tutti estesi (per esempio, ax diventa eax).

	*/

	asm {   push bp

		db 0x66; mov dx, word ptr b   // passa in edx il delta_y (ovvero b).
		db 0x66; mov ax, word ptr x2  // passa in eax il limite superiore del ciclo (come coordinata x d'arrivo, ovvero x2).
		db 0x66; mov bp, word ptr a } // passa in ebp il delta_x (ovvero a).

_do:	asm {   mov bx, [di]
		shl bx, 2
		add bx, [di]
		shl bx, 4
		db 0x66; add word ptr global_y, dx
		mov cx, [si]
		shr cx, 3
		add bx, cx
		mov cl, [si]
		mov ch, 128
		db 0x66; add word ptr global_x, bp
		and cl, 7
		shr ch, cl
		or es:[bx+4], ch
		db 0x66; cmp word ptr global_x, ax
		jb _do }

	asm 	pop bp;
}

/* Tracciamento bastoncini (parte 3d). */

void stick3d (double p_x, double p_y, double p_z,
	      double x, double y, double z)
{
	long   fx, fy, lx, ly;

	double rx, ry, rz, z2, kk;

	double p_rx, p_ry, p_rz, p_z2, diff;

	p_x -= cam_x;
	p_y -= cam_y;
	p_z -= cam_z;

	p_z2 = p_z * cosbeta - p_x * sinbeta;
	p_rz = p_z2 * cosalfa + p_y * sinalfa;

	x -= cam_x;
	y -= cam_y;
	z -= cam_z;

	z2 = z * cosbeta - x * sinbeta;
	rz = z2 * cosalfa + y * sinalfa;

	if (rz<min_dist&&p_rz<min_dist) return;

	p_rx = p_x * cosbeta_x_ob + p_z * sinbeta_x_ob;
	p_ry = p_y * cosalfa_x_ob - p_z2 * sinalfa_x_ob;

	rx = x * cosbeta_x_ob + z * sinbeta_x_ob;
	ry = y * cosalfa_x_ob - z2 * sinalfa_x_ob;

	/* Conversione punti alle spalle dell'osservatore rispetto al piano
	   dello schermo. */

	if (rz<min_dist) {
		diff = rz-p_rz;
		if (diff<-mindiff || diff>mindiff) {
			kk = (min_dist-p_rz) / diff;
			rx = kk * (rx-p_rx) + p_rx;
			ry = kk * (ry-p_ry) + p_ry;
		}
		rz = min_dist;
	}
	else {
		if (p_rz<min_dist) {
			diff = p_rz-rz;
			if (diff<-mindiff || diff>mindiff) {
				kk = (min_dist-rz) / diff;
				p_rx = kk * (p_rx-rx) + rx;
				p_ry = kk * (p_ry-ry) + ry;
			}
			p_rz = min_dist;
		}
	}

	// Prospettiva.

	fx = rx / rz; fy = ry / rz;
	lx = p_rx / p_rz; ly = p_ry / p_rz;

	/* Ottimizza il tracciamento del "Segmento", tagliando via le parti
	   che di sicuro non si vedono. Questa parte non ä efficace se
	   i punti di partenza e di arrivo fanno sç che il segmento non
	   intersechi alcun lato dell'area visibile: ma in questo peculiare
	   caso le linee vengono escluse dai controlli successivi. */

	if (fx<stk_lbx) {
		diff = fx-lx;
		if (diff<-mindiff || diff>mindiff) {
			kk = (stk_lbx-lx) / diff;
			fy = kk * (fy-ly) + ly;
			fx = stk_lbx;
		}
	}
	if (lx<stk_lbx) {
		diff = lx-fx;
		if (diff<-mindiff || diff>mindiff) {
			kk = (stk_lbx-fx) / diff;
			ly = kk * (ly-fy) + fy;
			lx = stk_lbx;
		}
	}
	if (fy<stk_lby) {
		diff = fy-ly;
		if (diff<-mindiff || diff>mindiff) {
			kk = (stk_lby-ly) / diff;
			fx = kk * (fx-lx) + lx;
			fy = stk_lby;
		}
	}
	if (ly<stk_lby) {
		diff = ly-fy;
		if (diff<-mindiff || diff>mindiff) {
			kk = (stk_lby-fy) / diff;
			lx = kk * (lx-fx) + fx;
			ly = stk_lby;
		}
	}
	if (fx>stk_ubx) {
		diff = fx-lx;
		if (diff<-mindiff || diff>mindiff) {
			kk = (stk_ubx-lx) / diff;
			fy = kk * (fy-ly) + ly;
			fx = stk_ubx;
		}
	}
	if (lx>stk_ubx) {
		diff = lx-fx;
		if (diff<-mindiff || diff>mindiff) {
			kk = (stk_ubx-fx) / diff;
			ly = kk * (ly-fy) + fy;
			lx = stk_ubx;
		}
	}
	if (fy>stk_uby) {
		diff = fy-ly;
		if (diff<-mindiff || diff>mindiff) {
			kk = (stk_uby-ly) / diff;
			fx = kk * (fx-lx) + lx;
			fy = stk_uby;
		}
	}
	if (ly>stk_uby) {
		diff = ly-fy;
		if (diff<-mindiff || diff>mindiff) {
			kk = (stk_uby-fy) / diff;
			lx = kk * (lx-fx) + fx;
			ly = stk_uby;
		}
	}

	if (fx==lx&&fy==ly) return; // Esclude le linee costituite da un punto solo.

	if (fy<stk_lby||ly<stk_lby) return; // Esclude le linee che mai e poi mai si vedranno.
	if (fy>stk_uby||ly>stk_uby) return;
	if (fx<stk_lbx||lx<stk_lbx) return;
	if (fx>stk_ubx||lx>stk_ubx) return;

	Stick (fx+320, fy+240, lx+320, ly+240);
}

/* Tracciamento di un punto.
   La funzione, altamente specifica, assolutamente inesportabile,
   usa la pagina nascosta ADAPTED, ed alcune variabili devono essere
   precaricate: shift (=y*80) ed x2d, ovvero x.
   Oltretutto, ä necessario caricare shift con un valore pió grande del
   necessario di 4 unitÖ: perchÇ questo ä l'allineamento fisso di adapted.
   Ma c'ä ancora qualcosa: in es si deve precaricare il valore seg_adapted. */

void plot ()
{
	asm {	mov bx, x2d // istruz. sfalsate per motivi di accoppiabilitÖ
		mov cx, x2d
		shr bx, 3
		and cl, 7
		mov dl, 128
		mov di, shift
		shr dl, cl
		or es:[bx+di], dl }
}

/* Scrive una stringa sul video, in caratteri di piccola dimensione ma
   altamente leggibili perchÇ intensi ed essenziali: ne vado fiero.
   Oltretutto, pur usando un bitmapped font, ä discretamente rapida
   (probabilmente molto pió veloce delle funzioni rom) perchÇ il suo
   nucleo ä stato ottimizzato pió volte.
   Non credo le rimanga molto da ottimizzare. */

void Testo (int x, int y, char *messaggio, unsigned stile)
{
	// Funzione modificata per agire con VGA hi-res su pagina nascosta.

	unsigned r;

	if (glow) return;

	if (stile==ALL_DX|stile==CENTRATO) {
		r = 0;
		_CX = 0;
		while (messaggio[_CX]) {
			if (messaggio[_CX]!='|')
				r += 6;
			else
				r ++;
			_CX++;
		}
		if (stile==ALL_DX)
			x += 638 - r;
		else
			x += 320 - (r >> 1);
	}

	asm {	mov si, word ptr messaggio
		les di, dword ptr adapted
		mov ax, y
		shl ax, 2
		add ax, y
		shl ax, 4
		add di, ax
		mov al, [si] }
_while: asm {   cmp al, 0
		jne next }

		return;

next:	asm {   push si
		xor ah, ah
		sub ax, 32
		mov bx, x
		shr bx, 3
		mov cx, x
		and cl, 7
		mov dl, 6
		mov si, ax
		add si, ax
		add si, ax
		shl si, 1 }
lineup: asm {	mov ah, byte ptr font[si]
		xor al, al
		shr ax, cl
		xchg ah, al
		or es:[bx+di], ax
		add bx, 80
		inc si
		dec dl
		jnz lineup
		pop si
		inc si
		cmp byte ptr [si], '|'
		jne forw
		inc x
		mov al, [si-1]
		jmp _while }
forw:	asm {	add x, 6
		mov al, [si]
		jmp _while }
}

/* Alloca la memoria necessaria. */

char m_init ()
{
	if (m) farfree (m);
	m = (float far *)  farmalloc ((long)ms*(long)M_GAL*4);
	if (!m) return (0);

	if (x) farfree (x);
	x = (double far *)  farmalloc ((long)ms*(long)M_GAL*8);
	if (!x) return (0);
	if (y) farfree (y);
	y = (double far *)  farmalloc ((long)ms*(long)M_GAL*8);
	if (!y) return (0);
	if (z) farfree (z);
	z = (double far *)  farmalloc ((long)ms*(long)M_GAL*8);
	if (!z) return (0);

	if (vx) farfree (vx);
	vx = (double far *) farmalloc ((long)ms*(long)M_GAL*8);
	if (!vx) return (0);
	if (vy) farfree (vy);
	vy = (double far *) farmalloc ((long)ms*(long)M_GAL*8);
	if (!vy) return (0);
	if (vz) farfree (vz);
	vz = (double far *) farmalloc ((long)ms*(long)M_GAL*8);
	if (!vz) return (0);

	if (adapted) farfree (adapted);
	adapted = (unsigned char far *) farmalloc (38400);
	if (!adapted) return (0);

	asm {
		push es
		push ax
		les ax, dword ptr m
		mov seg_m, es
		les ax, dword ptr x
		mov seg_x, es
		les ax, dword ptr y
		mov seg_y, es
		les ax, dword ptr z
		mov seg_z, es
		les ax, dword ptr vx
		mov seg_vx, es
		les ax, dword ptr vy
		mov seg_vy, es
		les ax, dword ptr vz
		mov seg_vz, es
		les ax, dword ptr adapted
		mov seg_adapted, es
		pop ax
		pop es
	}

	return (1);
}

/* Calcola le coordinate 2d di un punto in 3d. */

char conv (double x, double y, double z)
{
	x -= cam_x;
	y -= cam_y;
	z -= cam_z;

	rx = x * cos(beta) + z * sin(beta);
	z2 = z * cos(beta) - x * sin(beta);

	rz = z2 * cos(alfa) + y * sin(alfa);

	if (rz<1) return(0); // Il punto non si vede. Inutile continuare.

	ry = y * cos(alfa) - z2 * sin(alfa);

	w = ob / rz;

	x2d = rx * w + 320;
	y2d = ry * w + 240;

	return(1);
}

/* Orienta la videocamera verso il punto (see_x; see_y; see_z).
   Come implementazione ä penosa: se c'ä una cosa che potrï migliorare,
   forse, sarÖ questa funzione. Dunque: non ho trovata traccia di alcuna
   funzione matematica che permetta di trovare l'angolazione giusta in base
   alle coordinate e perciï ho dovuto cercare il punto giusto facendo
   scandire l'angolo giro per entrambi gli angoli.
   Ma i risultati sono a dir poco deludenti: innanzi tutto la funzione ä
   lenta (dovendo scandire tutti i 360¯ per due angoli, anche se in modo
   frattale, non puï essere altrimenti), e poi non inquadra sempre in modo
   perfetto. Sperimentalmente, ho visto che la si deve chiamare almeno
   due volte di seguito per esser certi che il punto venga inquadrato.

   Ho trovato la funzione atan (arcotangente), ma non so se riscrivere
   watch() in modo che la usi.

   Un po' perchÇ non ne ho voglia, e poi ä una funzione limitata, che
   potrebbe dare dei problemi di dominio. */

void watch ()
{
	rx     = cam_x - see_x;
	ry     = cam_y - see_y;
	rz     = cam_z - see_z;
	w      = sqrt (rx*rx+ry*ry+rz*rz);

	see_x -= cam_x;
	see_y -= cam_y;
	see_z -= cam_z;

	f1     = 0;
	f2     = 0;

	f3     = 2 * M_PI;

	while (f3 > 0.01) {

		q       = 1e+111;
		ww      = f3 / 10;

		sinbeta = sin(beta);
		cosbeta = cos(beta);

		for (alfa=f2-f3; alfa<f2+f3; alfa+=ww) {
			rx  = - w * sinbeta * cos(alfa) - see_x;
			ry  = + w * sin(alfa) 	        - see_y;
			rz  = + w * cosbeta * cos(alfa) - see_z;
			k   = rx*rx + ry*ry + rz*rz;
			if (k<q) {
				q = k;
				f2 = alfa;
			}
		}

		alfa    = f2;

		ry      = w * sin(alfa) - see_y;
		ry     *= ry;

		cosalfa = cos(alfa);

		q       = 1e+111;
		ww      = f3 / 10;

		for (beta=f1-f3; beta<f1+f3; beta+=ww) {
			rx  = - w * sin(beta) * cosalfa - see_x;
			rz  = + w * cos(beta) * cosalfa - see_z;
			k   = rx*rx + ry + rz*rz;
			if (k<q) {
				q = k;
				f1 = beta;
			}
		}

		beta = f1;

		f3 /= 5;

	}

	q = sqrt (see_x*see_x + see_y*see_y + see_z*see_z);
	conv (see_x + cam_x, see_y + cam_y + q, see_z + cam_z);
	if (y2d<240) {
		alfa += M_PI;
		beta += M_PI;
	}

	see_x += cam_x;
	see_y += cam_y;
	see_z += cam_z;
}

/* Memorizza la posizione e l'angolo di osservazione. */

void mem ()
{
	mcam_x = cam_x;
	mcam_y = cam_y;
	mcam_z = cam_z;

	m_alfa = alfa;
	m_beta = beta;
}

/* Ripristina i dati salvati dalla f. precedente. */

void rcm ()
{
	cam_x = mcam_x;
	cam_y = mcam_y;
	cam_z = mcam_z;

	alfa = m_alfa;
	beta = m_beta;
}

/* Fa avanzare la video camera nella direzione di volo. */

void forward (double delta)
{
	cam_x -= delta * sin(beta) * cos(alfa);
	cam_z += delta * cos(beta) * cos(alfa);
	cam_y += delta * sin(alfa);
}

/* Carica le aree dati riguardanti le stelle. */

void s_init (unsigned pan, char speeds)
{
	long double D, v, mA;
	long double px, py, pz;
	long double s1, c1, s2, c2;
	long double ang_1, ang_2, zeta;

	srand (X[pan]);

	i = pan * ms;
	D = d[pan] * A;
	mA = s[pan] * M * G;

	for (aux=0; aux<ms; i++, aux++) {

		m[i]  = (float)random(30000) * A;
		zeta  = D * (double)random(10000) / 20000;
		zeta -= D * (double)random(10000) / 20000;
		if (zeta<0) zeta = -zeta;
		zeta += A;

		switch (t[pan]) {
			case Ex:
				ang_1 = M_PI * (double)random(10000) / 10000;
				ang_2 = M_PI * (double)random(20000) / 10000;
				break;
			case Sa:
				ang_1  = M_PI * (double)random(15000) / 1E6;
				ang_1 -= M_PI * (double)random(15000) / 1E6;
				ang_2  = M_PI * (double)random(30000) / 1E6;
				ang_2 -= M_PI * (double)random(30000) / 1E6;
				if (zeta > D / 15) {
					ang_1 *= D / zeta;
					ang_2 *= D / zeta;
				}
				else {
					ang_1 *= 200;
					ang_2 *= 100;
				}
				if (random(2)) ang_2 += M_PI;
				ang_2 += M_PI_2;
				break;
			case Sc:
				ang_1  = M_PI * (double)random(15000) / 5E6;
				ang_1 -= M_PI * (double)random(15000) / 5E6;
				ang_2  = M_PI * (double)random(30000) / 5E6;
				ang_2 -= M_PI * (double)random(30000) / 5E6;
				if (zeta > D / 15) {
					ang_1 *= D / zeta;
					ang_2 *= D / zeta;
				}
				else {
					ang_1 *= 200;
					ang_2 *= 100;
				}
				if (random(2)) ang_2 += M_PI;
				ang_2 += M_PI_2;
				break;
			case S4:
				ang_1  = M_PI * (double)random(15000) / 2E6;
				ang_1 -= M_PI * (double)random(15000) / 2E6;
				ang_2  = M_PI * (double)random(30000) / 2E6;
				ang_2 -= M_PI * (double)random(30000) / 2E6;
				if (zeta > D / 15) {
					ang_1 *= D / zeta;
					ang_2 *= D / zeta;
				}
				else {
					ang_1 *= 200;
					ang_2 *= 100;
				}
				if (random(2)) ang_2 += M_PI;
				if (random(2)) ang_2 += M_PI_2;
		}

		ang_1 += incl[pan];
		ang_2 += tilt[pan];

		s1 = sin(ang_1);
		s2 = sin(ang_2);
		c1 = cos(ang_1);
		c2 = cos(ang_2);

		x[i]  = zeta * s2;
		y[i]  = - zeta * c2 * s1;
		z[i]  = zeta * c2 * c1;

		if (speeds) {
			v     = sqrtl (mA / zeta) * single_step;
			w     = zeta * c2 - v * s2;
			px    = v * c2 + zeta * s2;
			py    = - w * s1;
			pz    = w * c1;
			if (acws[pan]) {
				vx[i] = x[i] - px;
				vy[i] = y[i] - py;
				vz[i] = z[i] - pz;
			}
			else {
				vx[i] = px - x[i];
				vy[i] = py - y[i];
				vz[i] = pz - z[i];
			}
			vx[i] += gvx[pan];
			vy[i] += gvy[pan];
			vz[i] += gvz[pan];
		}

		x[i] += gx[pan];
		y[i] += gy[pan];
		z[i] += gz[pan];

	}
}

void g_init ()
{
	for (pan=0; pan<g; pan++)
		s_init (pan, 1);

	if (started&&!starting) {
		started = 0;
		follow = 0;
	}

	starting = 0;
}

/* Forze di gravitazione. */

void g_forces ()
{
	double ldrx, ldry, ldrz, division_limit = 1;
	double w = single_step * single_step * M * G;

	for (pan=0; pan<g; pan++) {
		_SI = ((pan * ms) << 3) + 4;
		asm {   lea bx, n
			mov di, pan
			shl di, 1
			mov cx, [bx+di] }
	glx_lp: asm {   mov dx, 0 }
	int_lp: asm {	mov di, dx
			shl di, 3
			mov es, seg_x
			fld qword ptr es:[si]
			fsub qword ptr gx[di]
			fst ldrx
			fmul ldrx
			mov es, seg_y
			fld qword ptr es:[si]
			fsub qword ptr gy[di]
			fst ldry
			fmul ldry
			mov es, seg_z
			fld qword ptr es:[si]
			fsub qword ptr gz[di]
			fst ldrz
			fmul ldrz
			faddp
			faddp
			fadd division_limit
			fstp q
			fld w
			fmul qword ptr s[di]
			fdiv q
			fld q
			fsqrt
			fdivp
			fstp k
			mov es, seg_vx
			fld qword ptr es:[si]
			fld ldrx
			fmul k
			fsubp st(1), st(0)
			fstp qword ptr es:[si]
			mov es, seg_vy
			fld qword ptr es:[si]
			fld ldry
			fmul k
			fsubp st(1), st(0)
			fstp qword ptr es:[si]
			mov es, seg_vz
			fld qword ptr es:[si]
			fld ldrz
			fmul k
			fsubp st(1), st(0)
			fstp qword ptr es:[si]
			inc dx
			cmp dx, g
			jnb nxt_g
			jmp int_lp }
	nxt_g:	asm {	add si, 8
			dec cx
			jz stars_end
			jmp glx_lp }
stars_end: }

	if (g>1) { // se vi ä pió d'una galassia...
	  for (aux=0; aux<g; aux++) { // ...tutte quante interagiscono...
		  for (pan=0; pan<g; pan++) { // ...con tutte le altre...
			if (aux!=pan) { // ...(tranne che con se stesse)
				ldrx = gx[aux] - gx[pan]; // componente x distanza
				ldry = gy[aux] - gy[pan]; // componente y distanza
				ldrz = gz[aux] - gz[pan]; // componente z distanza
				q    = ldrx*ldrx
				     + ldry*ldry
				     + ldrz*ldrz + A;	  // Distanza^2 (in metri).
				k = w * s[pan] / q;       // Accelerazione in m / step^2
				k /= sqrtl (q);
				gvx[aux] -= ldrx * k;	  // Incremento componente x veloc.
				gvy[aux] -= ldry * k;     //      -          -     y.
				gvz[aux] -= ldrz * k;     //      -          -     z.
			}
		  }
	  }
	}
}

/* Movimenti stellari. */

void tripper ()
{
	for (pan=0; pan<g; pan++) {
		_SI = ((pan * ms) << 3) + 4;
		asm {   lea bx, n
			mov di, pan
			shl di, 1
			mov cx, [bx+di] }
	pan_lp: asm {	mov es, seg_vx
			fld qword ptr es:[si]
			mov es, seg_x
			fadd qword ptr es:[si]
			fstp qword ptr es:[si]
			mov es, seg_vy
			fld qword ptr es:[si]
			mov es, seg_y
			fadd qword ptr es:[si]
			fstp qword ptr es:[si]
			mov es, seg_vz
			fld qword ptr es:[si]
			mov es, seg_z
			fadd qword ptr es:[si]
			fstp qword ptr es:[si]
			add si, 8
			dec cx
			jnz pan_lp }
		gx[pan] += gvx[pan];
		gy[pan] += gvy[pan];
		gz[pan] += gvz[pan];
	}

	TT += single_step;
}

/* Renderizza la scena sul video. */

void render ()
{
	unsigned count, index;

	for (pan=0; pan<g; pan++) {

		index = (pan * ms) << 3;

		asm {   lea bx, n
			mov di, pan
			shl di, 1
			mov ax, [bx+di]
			mov count, ax }

	pan_lp: asm     mov si, index

		// Traslazione dello spazio rispetto all'osservatore.

		asm {	mov es, seg_y
			fld qword ptr es:[si+4]
			fsub cam_y
			fstp _y
			mov es, seg_x
			fld qword ptr es:[si+4]
			fsub cam_x
			fst _x
			fmul sinbeta // Rotazione su asse Y.
			mov es, seg_z
			fld qword ptr es:[si+4]
			fsub cam_z
			fst _z
			fmul cosbeta
			fsubrp
			fst z2
			fmul cosalfa // Rotazione su asse X.
			fld _y
			fmul sinalfa
			faddp        // sullo stack c'ä ora rz.
			fcom min_dist
			fstsw ax
			sahf
			ja after
			ffree st(0)
			jmp useless } // Il punto non si vede (rz < min_dist): inutile continuare.

	after:  asm {   fld one // fld1 ä pió lenta (caso a 64 bit) di un ciclo.
			fdivrp
			fstp rz } // rz = 1 / rz: si moltiplica per dividere: + veloce.

		//	ry = (y * cosalfa_x_ob - z2 * sinalfa_x_ob) / rz;

		asm {	fld _y
			fmul cosalfa_x_ob
			fld z2
			fmul sinalfa_x_ob
			fsubp st(1), st(0)
			fmul rz
			fstp ry }

		_AX = (int)ry + 240;

		asm {	cmp ax, 3
			jnb ynbelow
			jmp useless }
	ynbelow:asm {	cmp ax, 476
			jb ybelow
			jmp useless }

		//	rx = (x * cosbeta_x_ob + z * sinbeta_x_ob) / rz;

	ybelow:	asm {   mov y2d, ax
			fld _x
			fmul cosbeta_x_ob
			fld _z
			fmul sinbeta_x_ob
			faddp
			fmul rz
			fstp rx }

		_AX = (int)rx + 320;

		asm {	cmp ax, 3
			jnb xnbelow
			jmp useless }
	xnbelow:asm {	cmp ax, 636
			jb xbelow
			jmp useless }

		// Tracciamento (differenziato, a seconda di distanza e massa).

		// z2d = m * ob / rz; (PoichÇ le masse sono premoltiplicate, m / rz).

	xbelow:	asm {   mov x2d, ax
			mov es, seg_m
			shr si, 1
			fld dword ptr es:[si+4]
			fmul rz
			fst z2d
			fcomp size_a
			fstsw ax
			sahf
			jb min_a
			jmp maj_a }

	min_a:	asm {	les di, adapted
			mov ax, y2d
			shl ax, 2
			add ax, y2d
			shl ax, 4
			add di, ax
			mov bx, x2d
			shr bx, 3
			mov cx, x2d
			and cl, 7
			mov dl, 128
			shr dl, cl
			or es:[bx+di], dl }
		asm 	jmp useless

	maj_a:	asm {	fld z2d
			fcomp size_b
			fstsw ax
			sahf
			jb min_b
			jmp maj_b }

	min_b:	asm { 	mov ax, y2d
			shl ax, 2
			add ax, y2d
			shl ax, 4
			add ax, 4
			mov shift, ax
			mov es, seg_adapted }

			plot ();
		asm	dec x2d
			plot ();
		asm	add x2d, 2
			plot ();
		asm	sub shift, 80
		asm	dec x2d
			plot ();
		asm	add shift, 160
			plot ();
		asm	jmp useless

	maj_b:	asm { 	mov ax, y2d
			shl ax, 2
			add ax, y2d
			shl ax, 4
			add ax, 4
			mov shift, ax
			mov es, seg_adapted }

		asm {	fld z2d
			fcomp size_c
			fstsw ax
			sahf
			jb min_c
			jmp maj_c }

	min_c:		plot ();
		asm	sub x2d, 2
			plot ();
		asm	inc x2d
			plot ();
		asm	add x2d, 2
			plot ();
		asm	inc x2d
			plot ();
		asm	sub shift, 160
		asm	sub x2d, 2
			plot ();
		asm     add shift, 80
			plot ();
		asm	add shift, 160
			plot ();
		asm	add shift, 80
			plot ();
		asm	jmp useless

	maj_c:		plot ();
		asm	sub x2d, 3
			plot ();
		asm 	inc x2d
			plot ();
		asm	inc x2d
			plot ();
		asm	add x2d, 2
			plot ();
		asm	inc x2d
			plot ();
		asm	inc x2d
			plot ();
		asm	sub shift, 240
		asm	sub x2d, 3
			plot ();
		asm	add shift, 80
			plot ();
		asm	add shift, 80
		asm	dec x2d
			plot ();
		asm	inc x2d
			plot ();
		asm	inc x2d
			plot ();
		asm	add shift, 160
			plot ();
		asm	dec x2d
			plot ();
		asm	dec x2d
			plot ();
		asm     add shift, 80
		asm	inc x2d
			plot ();
		asm	add shift, 80
			plot ();

	useless: asm {	add index, 8
			dec count
			jz tonext
			jmp pan_lp }

	tonext: }
}

/* Simile a render, cambia molto l'uso dello schermo. */

void glowin_render ()
{
	int      pict_w;
	unsigned r, sw, count, index;

	for (pan=0; pan<g; pan++) {

		index = (pan * ms) << 3;

		asm {   lea bx, n
			mov di, pan
			shl di, 1
			mov ax, [bx+di]
			mov count, ax }

	pan_lp: asm     mov si, index

		// Traslazione dello spazio rispetto all'osservatore.

		asm {	mov es, seg_y
			fld qword ptr es:[si+4]
			fsub cam_y
			fstp _y
			mov es, seg_x
			fld qword ptr es:[si+4]
			fsub cam_x
			fst _x
			fmul sinbeta // Rotazione su asse Y.
			mov es, seg_z
			fld qword ptr es:[si+4]
			fsub cam_z
			fst _z
			fmul cosbeta
			fsubrp
			fst z2
			fmul cosalfa // Rotazione su asse X.
			fld _y
			fmul sinalfa
			faddp        // sullo stack c'ä ora rz.
			fcom min_dist
			fstsw ax
			sahf
			ja after
			ffree st(0)
			jmp useless } // Il punto non si vede (rz < min_dist): inutile continuare.

	after:  asm {   fld one // fld1 ä pió lenta (caso a 64 bit) di un ciclo.
			fdivrp
			fstp rz } // rz = 1 / rz: si moltiplica per dividere: + veloce.

		//	ry = (y * cosalfa_x_ob - z2 * sinalfa_x_ob) / rz;

		asm {	fld _y
			fmul cosalfa_x_ob
			fld z2
			fmul sinalfa_x_ob
			fsubp st(1), st(0)
			fmul rz
			fstp ry }

		_AX = (int)ry + 240;

		asm {	cmp ax, 21
			jnb ynbelow
			jmp useless }
	ynbelow:asm {	cmp ax, 458
			jb ybelow
			jmp useless }

		//	rx = (x * cosbeta_x_ob + z * sinbeta_x_ob) / rz;

	ybelow:	asm {   mov spy, ax
			fld _x
			fmul cosbeta_x_ob
			fld _z
			fmul sinbeta_x_ob
			faddp
			fmul rz
			fstp rx }

		_AX = (int)rx + 320;

		asm {	cmp ax, 21
			jnb xnbelow
			jmp useless }
	xnbelow:asm {	cmp ax, 618
			jb xbelow
			jmp useless }

	xbelow:	asm {   mov spx, ax
			mov es, seg_m
			shr si, 1
			fld dword ptr es:[si+4]
			fmul rz
			fcom size_g0
			fstsw ax
			sahf
			jnb maj_g0
			ffree st(0)
			jmp useless }
	maj_g0:	asm {	fcom size_g1
			fstsw ax
			sahf
			jb min_g1
			ffree st(0)
			fld size_g1 }

	/* Traccia la stella in modo glowin', ovvero come un piccolo
	   diamante luminoso la cui sagoma ä contenuta nel buffer "gws". */

	min_g1: asm {   fistp pict_w
			mov si, pict_w
			shl si, 1
			db 0x8b, 0x84
			dw offset gws_szk
			//mov ax, word ptr gws_szk[si]
			mov sw, ax
			shr ax, 1
			sub spx, ax
			sub spy, ax
			mov cx, p_align
			mov ax, spy
			xor dx, dx
			div cx
			mov cr_sp_page, ax
			shl ax, 8
			mov r, ax
			mov ax, dx
			mov dx, word ptr screenlength
			mul dx
			add ax, spx
			mov di, ax
			cmp ax, r
			jnb if1
			dec cr_sp_page }
	   if1: asm {	sub di, r
			db 0x8b, 0xb4
			dw offset gws_offset
			//mov si, word ptr gws_offset[si]
			mov bx, seg_desperation
			mov dx, cr_sp_page
			shl dx, 12
			add bx, dx
			mov es, bx
			mov dx, 640
			sub dx, sw
			mov ch, byte ptr sw }
	  line:	asm 	mov ah, byte ptr sw
	 pixel:	asm {   cmp byte ptr gws[si], 0
			je nop
			mov al, byte ptr gws[si]
			add es:[di], al
			jnc nop
			mov byte ptr es:[di], 0xff }
	   nop:	asm {   inc si
			add di, 1
			jnc smpage
			add bx, 0x1000
			mov es, bx }
	smpage:	asm {	dec ah
			jnz pixel }
		asm {   add di, dx
			jnc right
			add bx, 0x1000
			mov es, bx }
	 right:	asm {	dec ch
			jnz line }

       useless: asm {	add index, 8
			dec count
			jz tonext
			jmp pan_lp }

	tonext: }
}

/* Funzione per ridefinire i colori.
   Produce una tavola di 16 grigi. */

void grays (char h)
{
	char index[15] = { 1, 2, 3, 4, 5, 0x14, 7, 0x38, 0x39,
			   0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f };

	asm {
		pusha
		mov cx, 15
		xor bh, bh
		mov dl, 7
	}
ridef:  asm {
		push cx
		mov di, 15
		sub di, cx
		mov bl, byte ptr index[di]
		push dx
		xor ah, ah
		xor dh, dh
		mov al, dl
		mov dl, h
		mul dx
		pop dx
		shr ax, 6
		mov dh, al
		mov ch, al
		mov cl, al
		mov ax, 1010h
		int 10h
		add dl, 4
		pop cx
		loop ridef
		popa
	}
}

/* Disegna una stellina per la modalitÖ CCD. */

void ccd_star ()
{
	int a, b, c, cr_exp;

	unsigned mag, x, y;

	float f;

	asm {
		fld _x
		fld cam_x
		fsubp st(1), st(0)
		fstp _x
		fld _y
		fld cam_y
		fsubp st(1), st(0)
		fstp _y
		fld _z
		fld cam_z
		fsubp st(1), st(0)
		fst _z
		fld cosbeta
		fmulp
		fld _x
		fld sinbeta
		fmulp
		fsubp st(1), st(0)
		fst z2
		fld cosalfa
		fmulp
		fld _y
		fld sinalfa
		fmulp
		faddp
		fld min_dist
		fcomp
		fstsw ax
		sahf
		jbe after
		ffree st(0)
	}

	return;

after:  asm {
		fld one
		fdivrp
		fstp rz
		fld _y
		fld cosalfa_x_ob
		fmulp
		fld z2
		fld sinalfa_x_ob
		fmulp
		fsubp st(1), st(0)
		fld rz
		fmulp
		fstp ry
		fld _m
		fld rz
		fmulp
		fstp z2d
	}

	y2d = (int)ry + 240;
	mag = (int)(z2d * 10);
	if (mag>14) mag = 14;

	mag += mag_factor;

	if (!mag)	   return;
	if (y2d<mag)       return;
	if (y2d>479 - mag) return;

	asm {
		fld _x
		fld cosbeta_x_ob
		fmulp
		fld _z
		fld sinbeta_x_ob
		fmulp
		faddp
		fld rz
		fmulp
		fstp rx
	}

	x2d = (int)rx + 320;

	if (x2d<mag)       return;
	if (x2d>639 - mag) return;

	b = -mag;
	y = y2d - mag;
	while (y <= y2d + mag) {
		a = -mag;
		x = x2d - mag;
		while (x <= x2d + mag) {
			f = sqrt (a*a + b*b);
			if (f<mag) {
				c = exposure - (f * exposure / mag);
				if (c>1) {
					asm {
						mov ah, 13
						mov cx, x
						mov dx, y
						xor bh, bh
						int 0x10
						xor ah, ah
						add c, ax
					}
					if (c>14) c = 14;
					asm {
						mov al, byte ptr c
						mov ah, 12
						mov cx, x
						mov dx, y
						xor bh, bh
						int 0x10
					}
				}
			}
			x++;
			a++;
		}
		y++;
		b++;
	}

	if (_m < crlim) return;
	if (mag<7) return;

	mag *= 2;
	cr_exp = exposure * 2;

	if (y2d<mag)       return;
	if (y2d>479 - mag) return;
	if (x2d<mag)       return;
	if (x2d>639 - mag) return;

	// Disegna la diffrazione della crociera (simulata).

	b = -mag;
	y = y2d - mag;
	while (y <= y2d + mag) {
		if (b>0)
			c = cr_exp - (b * cr_exp / mag);
		else
			c = cr_exp - ( -b * cr_exp / mag);
		if (c>1) {
			asm {
				mov ah, 13
				mov cx, x2d
				mov dx, y
				xor bh, bh
				int 0x10
				xor ah, ah
				add c, ax
			}
			if (c>14) c = 14;
			asm {
				mov al, byte ptr c
				mov ah, 12
				mov cx, x2d
				mov dx, y
				xor bh, bh
				int 0x10
			}
		}
		y++;
		b++;
	}

	a = -mag;
	x = x2d - mag;
	while (x <= x2d + mag) {
		if (a>0)
			c = cr_exp - (a * cr_exp / mag);
		else
			c = cr_exp - ( -a * cr_exp / mag);
		if (c>1) {
			asm {
				mov ah, 13
				mov cx, x
				mov dx, y2d
				xor bh, bh
				int 0x10
				xor ah, ah
				add c, ax
			}
			if (c>14) c = 14;
			asm {
				mov al, byte ptr c
				mov ah, 12
				mov cx, x
				mov dx, y2d
				xor bh, bh
				int 0x10
			}
		}
		x++;
		a++;
	}
}

/* Produce uno schermo simile a quello che si potrebbe ottenere con un CCD.
   I parametri per la funzione di tracciamento stelle sono:
	mag_factor ... sfocatura immagine (0 = a fuoco): per default ä 3;
	exposure ..... tempo esposizione in minuti: solitamente da 3 a 6. */

void ccd_image ()
{
	pclear (adaptor, 0);
	grays (64);

	unsigned count, index;

	for (pan=0; pan<g; pan++) {
		crlim = A * (100 - (double)pcrlim) * 300;
		index = (pan * ms) << 2;
		asm {   lea bx, n
			mov di, pan
			shl di, 1
			mov ax, [bx+di]
			mov count, ax }
	pan_lp: asm {   mov si, index
			mov es, seg_m
			fld dword ptr es:[si+4]
			fstp _m
			shl si, 1
			mov es, seg_x
			fld qword ptr es:[si+4]
			fstp _x
			mov es, seg_y
			fld qword ptr es:[si+4]
			fstp _y
			mov es, seg_z
			fld qword ptr es:[si+4]
			fstp _z }
			ccd_star ();
		asm {	add index, 4
			dec count
			jnz pan_lp }
	}
}

/* Gestione bliper. */

unsigned cik[20] = { 1111, 11000, 15555, 12000, 14000, 1333, 1100, 180, 140, 120,
		     1111, 11000, 15555, 12000, 14000, 1333, 1100, 180, 140, 120 };

void blip (unsigned *seq)
{
	unsigned l = sizeof (seq) * 2;

	for (_SI = 0; _SI < l; _SI += 2) {
	asm {   push dx
		push ax
		push bx
		mov bx, word ptr seq[si]
		mov dx, 12h
		mov ax, 34ddh
		div bx
		mov bx, ax
		in al, 61h
		test al, 3
		jne neq
		or al, 3
		out 61h, al
		mov al, 0b6h
		out 43h, al }
neq:    asm {   mov al, bl
		out 42h, al
		mov al, bh
		out 42h, al
		pop bx
		pop ax
		pop dx }
	}

	asm {
		push ax
		in al, 61h
		and al, 0fch
		out 61h, al
		pop ax
	}
}

/* Tenta di ancorare l'osservatore alla stella inquadrata pió vicina
   e pió centrale sullo schermo. Dico "tenta" perchÇ non ä detto che ci siano
   per forza stelle inquadrate nel mirino. */

void join ()
{
	j = anchor;
	anchor = -1;
	ww = 111e+111;

	for (pan=0; pan<g; pan++)
		for (aux=0; aux<n[pan]; aux++) {
			i  = ms * pan + aux;
			if (conv (x[i], y[i], z[i])) {
				rx = x2d - 320;
				ry = y2d - 240;
				k  = sqrt (rx*rx + ry*ry);
				if (k<25) {
					if (i!=j) {
						k += rz;
						if (k<ww) {
							delta_pcam = 0;
							anchor = i;
							ww = k;
						}
					}
				}
			}
		}
}

/* Calcola le velocitÖ per le galassie destinate allo scontro,
   e ripristina le posizioni specificate per l'inizio della simulazione. */

void thrower ()
{
	double _camx, _camy, _camz;

	errorflag = -1;

	for (pan=0; pan<g; pan++) {
		gx[pan] = bgx[pan];
		gy[pan] = bgy[pan];
		gz[pan] = bgz[pan];
		if (trgt[pan] == pan) {
			gvx[pan] = 0;
			gvy[pan] = 0;
			gvz[pan] = 0;
		}
		else {
			_alfa = alfa;
			_beta = beta;
			_camx = cam_x;
			_camy = cam_y;
			_camz = cam_z;
			cam_x = gx[pan];
			cam_y = gy[pan];
			cam_z = gz[pan];
			see_x = bgx[trgt[pan]];
			see_y = bgy[trgt[pan]];
			see_z = bgz[trgt[pan]];
			watch ();
			watch ();
			alfa += ashf[pan];
			beta += bshf[pan];
			gvx[pan] = - vlct[pan] * sin(beta) * cos(alfa);
			gvz[pan] =   vlct[pan] * cos(beta) * cos(alfa);
			gvy[pan] =   vlct[pan] * sin(alfa);
			alfa  = _alfa;
			beta  = _beta;
			cam_x = _camx;
			cam_y = _camy;
			cam_z = _camz;
		}
	}
}

/* Produce lo schermo di presentazione di StarWays 1997. */

void intro ()
{
	int fh;
	char c;

	initvesa (_640);
	if (video_mode[0]!=_640) { // vesa non disponibile, allora...
		_320_200_256 ();   // ...mcga.
		tavola_colori (gpal, 0, 256, 0, 0, 0);
		fh = _open ("STARWAYS.FLR", 0);
		if (fh>-1) {
			_read (fh, gpal, 768);
			_read (fh, adaptor, 64000);
			_close (fh);
		}
	}
	else {
		pselect (0);
		tavola_colori (gpal, 0, 256, 0, 0, 0);
		asm {   les di, dword ptr adaptor
			mov cx, 16384
			xor si, si
			xor dx, dx
			db 0x66; xor ax, ax }
	clr_pg: asm {  	db 0xf3 // macro: rep stosd
			db 0x66
			db 0xab
			inc dx
			mov ax, 0x4f05
			mov bx, 1
			int 0x10
			mov cx, 16384
			db 0x66; xor ax, ax
			cmp dx, 4
			ja  clr_st
			jne clr_pg
			mov cx, 11263
			jmp clr_pg }
	clr_st: fh = _open ("STARWAYS.FLR", 0);
		if (fh>-1) {
			_read (fh, gpal, 768);
			_close (fh);
		}
		flr_on_vesa ("STARWAYS.FLR", 160, 140, 320, 200, 768);
	}

	c = 0;
	do {
		c++;
		tci = clock ();
		tavola_colori (gpal, 0, 256, c, c, c);
		while ((c%2)&&tci==clock());
	} while (c!=63);
	tci = clock ();
	while (clock()-tci<20);
	c = 63;
	do {
		c--;
		tci = clock ();
		tavola_colori (gpal, 0, 256, c, c, c);
		while ((c%2)&&tci==clock());
	} while (c);

	fh = _open ("GLOWSTAR.TAB", 0);
	if (fh>-1) {
		_read (fh, gpal, 768);
		_close (fh);
	}
}

void puntamento_nuclei_galattici ()
{
	msg = 3;

	k = 1e+308;
	for (pan=0; pan<g; pan++) {
		if (conv (gx[pan], gy[pan], gz[pan])) {
			rx = x2d - 320;
			ry = y2d - 240;
			ww = sqrt (rx*rx + ry*ry);
			if (ww<75) {
				rx = gx[pan] - cam_x;
				ry = gy[pan] - cam_y;
				rz = gz[pan] - cam_z;
				ww = rx*rx + ry*ry + rz*rz;
				if (ww<k) {
					k = ww;
					msg = 2;
					sglx = pan;
				}
			}
		}
	}

	if (glxpar==2) {
		sprintf (input, "%1.0f", vlct[sglx] / S);
		ic = strlen (input) - 1;
	}

	if (glxpar==5) {
		sprintf (input, "%d", n[sglx]);
		ic = strlen (input) - 1;
	}
}

/* Traccia una bolla testuale di lavoro. */

void bolla (int x_centro, unsigned y_centro, double raggio)
{
	if (glow) return;

	i = -1;

	ww = M_PI_4 / raggio;

	for (w = -M_PI_2; w < M_PI_2; w += ww) {
		shift = ((long)(raggio * sin(w)) + y_centro) * 80;
		if (shift!=i) {
			x2d = (x_centro - (int)(raggio * cos(w))) >> 3;
			_DX = (((int)(raggio * cos(w)) * 2) >> 3) + 2;
			asm {	les di, adapted
				add di, shift
				add di, x2d
				db 0x66
				xor ax, ax
				mov cx, dx
				shr cx, 2
				db 0xf3 // macro: rep stosd
				db 0x66
				db 0xab
				mov cx, dx
				and cx, 3
			    rep stosb }
			i = shift;
		}
	}

	k = bbangle + 2 * M_PI;
	ww = 2 * M_PI / raggio;

	asm mov es, seg_adapted;

	for (w = bbangle; w < k; w += ww) {
		shift = ((long)(raggio * sin(w)) + y_centro) * 80 + 4;
		x2d = x_centro + (int)(raggio * cos(w));
				       plot ();
			     x2d += 2; plot ();
		shift += 80; x2d --;   plot ();
		shift += 80; x2d --;   plot ();
			     x2d += 2; plot ();
	}
}

void store_input_value ()
{
	if (msg) {
		if (glxpar==2) vlct[sglx] = atof(input) * S;
		if (glxpar==4) {
			if (subglxpar==1) d[sglx] = atof(input);
			if (subglxpar==2) s[sglx] = atof(input) * 1E9;
		}
		if (glxpar==5) {
			n[sglx] = atoi (input);
			if (n[sglx]==0) n[sglx] = 1;
			if (n[sglx]>ms) n[sglx] = ms;
		}
	}

	if (ccdpar) {
		if (ccdpar==1) exposure = atoi (input);
		if (ccdpar==2) mag_factor = atof (input);
		if (ccdpar==3) pcrlim = atoi (input);
		ccdpar = 0;
	}
}

char file[13*16];
char ext[6];

#define LOAD	0 // per il flag "purpose".
#define SAVE	1
#define DELETE	2

char *file_selector (char *title, char *prev_filename, char *extension,
		     char purpose)
{
	int i = 0, n = 0, tf = 0, nf = -1;
	int done, c, y, p, a = 0, b, f;

	struct ffblk ffblk;

	char reload = 1;
	char modif = 0;
	char chr = 0;

	pcs = fpcs;

	strcpy (ext, "*.");
	strcat (ext, extension);

	tf = 0;
	done = findfirst(ext,&ffblk,0);
	while (!done) {
		if (!strcmp(prev_filename, ffblk.ff_name)) n = tf;
		done = findnext(&ffblk);
		tf++;
	}

	strcpy (prev_filename, "FILENAME");

rewhile:do {
		if (tasto_premuto()) {
			i = attendi_pressione_tasto();
			if (i==0) {
				i = attendi_pressione_tasto();
				switch (i) {
					case 71: if (n) reload = 1;
						 n = 0;
						 break;
					case 79: if (n<tf) reload = 1;
						 n = tf;
						 break;
					case 80: if (n%15==14) reload = 1;
						 if (n<tf) n++;
						 break;
					case 72: if (n&&n%15==0) reload = 1;
						 if (n>0) n--;
						 break;
					case 81: n = a + 15;
						 if (n>tf)
							n = tf;
						 else
							reload = 1;
						 break;
					case 73: n -= 15;
						 if (n<0)
							n += 15;
						 else
							reload = 1;
						 break;
				}
				i = 0;
			}
			else {
				if (n==tf) {
					if (i==8) {
						if (!modif) {
							file[13*nf] = 0;
							chr = 0;
						}
						else {
							if (chr>0) {
								chr--;
								file[13*nf+chr] = 0;
							}
						}
						strcpy (prev_filename, file+13*nf);
					}
					else {
						if (i>32&&i<127&&chr<8) {
							if (i=='[') goto deny;
							if (i==']') goto deny;
							if (i=='<') goto deny;
							if (i=='>') goto deny;
							if (i=='.') goto deny;
							if (i==',') goto deny;
							if (i==';') goto deny;
							if (i==':') goto deny;
							if (i=='?') goto deny;
							if (i=='*') goto deny;
							if (i=='+') goto deny;
							if (i=='|') goto deny;
							if (i=='/') goto deny;
							if (i=='\\') goto deny;
							if (i>='a'&&i<='z') i -= 32;
							file[13*nf+chr+1] = 0;
							file[13*nf+chr] = i;
							modif = 1;
							chr++;
						  deny: strcpy (prev_filename,
								file+13*nf);
						}
					}
				}
			}
		}
		//
		a = (n/15) * 15;
		b = a + 15;
		if (reload) {
			tf = 0; c = 0;
			done = findfirst(ext,&ffblk,0);
			while (!done) {
				if (tf>=a&&tf<b) {
					ffblk.ff_name[strlen(ffblk.ff_name)-4] = 0;
					strcpy (file+13*c, ffblk.ff_name);
					c++;
				}
				done = findnext(&ffblk);
				tf++;
			}
			reload = 0;
			if (c<15) {
				strcpy (file+13*c, prev_filename);
				nf = c;
			}
			else
				nf = -1;
		}
		//
		pclear (adapted, 0);
		Testo (0,  0, "- F|I|L|E S|E|L|E|C|T|I|O|N  |B|U|B|B|L|E -", CENTRATO);
		Testo (0, 20, " |E|S|C| | @ CANCELS SELECTION AND RETURNS TO SIMULATION.", CENTRATO);
		Testo (0, 30, "E|N|T|E|R| @ CONFIRMS SELECTION AND EXECUTES.             ", CENTRATO);
		bolla (320, 240, 160);
		Testo (0, 134, title, CENTRATO);
		sprintf (tx, "LISTING FILES WITH EXTENSION: %s.", extension);
		Testo (0, 150, tx, CENTRATO);
		//
		p = 0;
		y = 180;
		for (c=a; c<b; c++) {
			if (c>tf) break;
			if (c==tf) {
				Testo (285, y, file+p, NORMALE);
				if (clock()%6<3) {
					if (c==n)
						Testo (285+6*strlen(file+p), y, "&", NORMALE);
					else
						Testo (284, y, file+p, NORMALE);
				}
				if (c==n) {
					if (purpose==SAVE)
						Testo (235, y, "N|E|W|>", NORMALE);
					else
						Testo (199, y, "S|P|E|C|I|F|Y|>", NORMALE);
				}
			}
			else {
				Testo (285, y, file+p, NORMALE);
				Testo (333, y, ".", NORMALE);
				Testo (339, y, extension, NORMALE);
				if (c==n) {
					Testo (250, y, ">|>|>", NORMALE);
					Testo (249, y, ">", NORMALE);
				}
			}
			p += 13;
			y += 8;
		}
		//
		if (!tf)
			strcpy (tx, "NO FILES OF THIS KIND.");
		else {
			if (tf>1)
				sprintf (tx, "%d FILES ACTUALLY STORED.", tf);
			else
				strcpy (tx, "1 FILE ACTUALLY STORED.");
		}
		Testo (0, 330, tx, CENTRATO);
		sprintf (tx, "P|A|G|E| %d O|F| %d.", n/15+1, tf/15+1);
		Testo (0, 346, tx, CENTRATO);
		//
		pcopy (adaptor, adapted);
		//
		bbangle += 0.0025;
		if (bbangle>M_PI_4) bbangle -= M_PI_4;
	} while (i!=27&&i!=13);

	if (i==27)
		return (NULL);
	else {
		i = -1;
		if (n==tf) {
			if (!file[13*nf])
				goto rewhile;
		}
		if (n%15) strcpy (file, file+13*(n%15));
		strcat (file, ".");
		strcat (file, extension);
		if (purpose==SAVE) {
			f = _open (file, 0);
			if (f>-1) {
				_close (f);
				Testo (0,   430, "A FILE WITH THAT NAME ALREADY EXISTS", CENTRATO);
				Testo (0,   440, "SHOULD IT BE O|V|E|R|W|R|I|T|T|E|N|?", CENTRATO);
				Testo (290, 460, "O|VERWRITE.", NORMALE);
				Testo (290, 470, "C|ANCEL.", NORMALE);
				pcopy (adaptor, adapted);
				i = attendi_pressione_tasto();
				while (tasto_premuto()) attendi_pressione_tasto();
				if (i!='o'&&i!='O') {
					reload = 1;
					goto rewhile;
				}
			}
		}
		if (purpose!=SAVE) {
			f = _open (file, 0);
			if (f==-1) {
				Testo (0, 445, ">>> F|I|L|E| N|O|T| A|C|C|E|S|S|I|B|L|E| <<<", CENTRATO);
				pcopy (adaptor, adapted);
				simwait ();
				reload = 1;
				goto rewhile;
			}
			_close (f);
		}
		if (purpose==DELETE) {
			f = _open (file, 0);
			if (f>-1) {
				_close (f);
				Testo (0, 430, "THE FILE WILL BE DELETED.", CENTRATO);
				Testo (260, 460, "D|ELETE THIS FILE.", NORMALE);
				Testo (260, 470, "C|ANCEL.", NORMALE);
				pcopy (adaptor, adapted);
				i = attendi_pressione_tasto();
				while (tasto_premuto()) attendi_pressione_tasto();
				if (i!='d'&&i!='D') {
					reload = 1;
					goto rewhile;
				}
			}
		}
		return (file);
	}
}

void last_exit () // Uscita d'emergenza.
{
	_80_25_C ();

	messaggio ("\r\nWarning:\r\n");
	messaggio ("\r\nThe file you tried to load is incomplete or damaged.");
	messaggio ("\r\nAt this point it would be advisable to verify the disk's structure");
	messaggio ("\r\nusing the CHKDSK or SCANDISK commands. This session of STARWAYS has");
	messaggio ("\r\nbeen stopped because its memory areas have been partly loaded with");
	messaggio ("\r\ninvalid data during the attempt to load that file.\r\n");
}

/* Sezione principale. */

void reach_your_dir ()
{
	int c;
	char d;

	c = lstri(_argv[0]) - 1;
	while (c>=0&&temptesto[c]!='\\') c--;

	if (c>=0) {
		if (temptesto[c-1] != ':')
			temptesto[c] = 0;
		else
			temptesto[c+1] = 0;
	}

	if (_argv[0][0]>='a'&&_argv[0][0]<='z')
		d = _argv[0][0] - 'a';
	else
		d = _argv[0][0] - 'A';

	asm {
		pusha
		mov dl, d
		mov ah, 0x0e
		int 0x21
		mov ah, 0x3b
		lea dx, temptesto
		int 0x21
		popa
	}
}

void main ()
{
	int fh, vptr;
	double vrif_y = 0;

	if (m_init ()) {
		thrower ();
		g_init ();
	}
	else {
		_80_25_C ();
		messaggio ("\r\nInsufficient conventional (DOS) memory.\r\n");
		return;
	}

	if (!test_and_init_mouse ()) {
		_80_25_C ();
		messaggio ("\r\nMouse or mouse driver not detected.\r\n");
		return;
	}

	harderr (handler);
	reach_your_dir ();

	// intro ();

	_640_480_16 (); // Apertura schermo grafico VGA.

	QUADWORDS = 9600; // Indicazione memoria video effettivamente usata.

	tci = clock ();

	while (1) { // Ciclo "a divinis". Non ridete.

		/* Campionamento mouse. */

		P = mpul;
		mpul = 0;
		mouse_input ();

		if (P==3&&mpul<3) {
			started = 0;
			starting = 0;
			if (msg==0) {
				if (glow)
					hylit = 1 - hylit;
				else
					join ();
			}
			else
				puntamento_nuclei_galattici ();
		}

		/* Campionamento tastiera. */

		key = -1;

		if (tasto_premuto()) {
			key = attendi_pressione_tasto();
			started = 0;
			starting = 0;
			if (exitflag) {
				if (!key) attendi_pressione_tasto();
				if (key=='N'||key=='n') exitflag = 0;
				if (key=='Y'||key=='y') {
					_80_25_C ();
					return;
				}
				key = -1;
			}
			if (glow) {
				if (!key) attendi_pressione_tasto();
				if (key!='v'&&key!='V'&&key!=27) key = -1;
				if (key==27) key = 'v';
			}
			if (key) {
				if (key==8) {
					input[ic] = 0;
					ic--;
					if (ic<0) {
						input[0] = '0';
						input[1] = 0;
						ic = 0;
					}
				}
				if (key>='A'&&key<='Z') key += 32;
				switch (key) {
					//
					case 'w':
						for (i=0; i<ms*M_GAL; i++) {
							vx[i] = -vx[i];
							vy[i] = -vy[i];
							vz[i] = -vz[i];
						}
						for (i=0; i<M_GAL; i++) {
							gvx[i] = -gvx[i];
							gvy[i] = -gvy[i];
							gvz[i] = -gvz[i];
						}
						break;
					//
					case 'a':
						if (!msg) {
							gas = 1 - gas;
							if (lens&&gas) lens = 0;
						}
						else {
							store_input_value ();
							subglxpar = 0;
							glxpar = 7;
						}
						break;
					case 'o':
						if (ccd_settings) {
							if (file_selector ("I|M|A|G|E F|I|L|E| |N|A|M|E|:", ccd_outfile, "BMP", SAVE)) {
								ccd_image ();
								pclear (adapted, 0);
								fh = _creat (file, 0);
								if (fh>-1) {
									_write (fh, bmpheader, 54);
									vptr = 0; while (vptr<256) {
										_write (fh, &vptr, 1);
										_write (fh, &vptr, 1);
										_write (fh, &vptr, 1);
										_write (fh, "\0", 1);
										vptr += 17;
									}
									vptr = 480; while (vptr) {
										asm {	les di, adapted
											mov dx, vptr
											dec dx
											mov si, 320
											xor cx, cx
											xor bh, bh }
									triga:	asm {	mov ah, 13
											int 0x10
											shl al, 4
											mov bl, al
											inc cx
											mov ah, 13
											int 0x10
											add al, bl
											stosb
											inc cx
											dec si
											jnz triga }
										_write (fh, adapted, 320);
										vptr--;
									}
									_close(fh);
									Testo (0, 470, "IMAGE STORED.", CENTRATO);
									strcpy (ccd_outfile, file);
								}
								else
									Testo (0, 470, "ERROR STORING FILE: DISK WRITE PROTECTED?", CENTRATO);
								QUADWORDS = 120;
								pcopy (adaptor+37600, adapted+37600);
								QUADWORDS = 9600;
								simwait ();
							}
						}
						else {
							lens = 1 - lens;
							if (lens&&gas) gas = 0;
						}
						break;
					case 'n':
						if (anchor!=65535) follow = +0;
						break;
					case 'x':
						if (ctrlkeys[0]&3) {
							alfa = 0;
							delta_alfa = 0;
							delta_beta = 0;
							cam_y = 0; cam_z = 0;
							if (cam_x>0)
								beta = M_PI_2;
							else
								beta = - M_PI_2;
						}
						else {
							if (anchor!=65535)
								follow = +1;
						}
						break;
					case 'y':
						if (ctrlkeys[0]&3) {
							beta = 0;
							delta_alfa = 0;
							delta_beta = 0;
							cam_x = 0; cam_z = 0;
							if (cam_y>0)
								alfa = - M_PI_2;
							else
								alfa = M_PI_2;
						}
						else {
							if (anchor!=65535)
								follow = -1;
						}
						break;
					case 'z':
						if (ctrlkeys[0]&3) {
							alfa = 0;
							delta_alfa = 0;
							delta_beta = 0;
							cam_x = 0; cam_y = 0;
							if (cam_z>0)
								beta = M_PI;
							else
								beta = 0;
						}
						else {
							if (anchor!=65535)
								follow = +2;
						}
						break;
					case 'f':
						if (ccd_settings) {
							store_input_value ();
							if (glxpar==2) glxpar = 0;
							if (subglxpar) subglxpar = 0;
							ccdpar = 2;
							sprintf (input, "%d", (int) mag_factor);
							ic = strlen (input) - 1;
							point = 0;
							sign = 0;
							ics = 1;
						}
						if (wait&&errorflag==-1) {
							tripper ();
							g_forces ();
							blip (cik);
							while (tasto_premuto())
								attendi_pressione_tasto();
						}
						break;
					case 'd':
						if (msg) {
							store_input_value ();
							glxpar = 4;
						}
						else
							anchor = -1;
						break;
					case 'g':
						grid = 1 - grid;
						break;
					case 'r':
						if (msg) {
							if (glxpar==1) acws[sglx] = 1 - acws[sglx];
							break;
						}
						TT = 0;
						thrower ();
						g_init ();
						break;
					case 's':
						wait = 1;
						break;
					case 'p':
						if (msg) break;
						wait = 0;
						break;
					case 'm':
						if (!msg) {
							TT = 0;
							anchor = -1;
							thrower ();
							g_init ();
							wait = 1;
							msg = 1;
						}
						else {
							thrower ();
							g_init ();
							msg = 0;
						}
						break;
					case '+':
						if (!msg) break;
						if (g<M_GAL) {
							msg     = 4;
							n[g]    = n[sglx];
							t[g]    = t[sglx];
							d[g]    = d[sglx];
							s[g]    = s[sglx];
							gx[g]   = cam_x;
							gy[g]   = cam_y;
							gz[g]   = cam_z;
							bgx[g]  = gx[g];
							bgy[g]  = gy[g];
							bgz[g]  = gz[g];
							trgt[g] = g;
							vlct[g] = 0;
							ashf[g] = 0;
							bshf[g] = 0;
							s_init (g, 0);
							sglx = g;
							g ++;
						}
						else
							msg = 5;
						break;
					case 't':
						if (msg) {
							store_input_value ();
							subglxpar = 0;
							glxpar = 3;
							cam_x = bgx[sglx];
							cam_y = bgy[sglx];
							cam_z = bgz[sglx];
							if (trgt[sglx]!=sglx) {
								see_x = bgx[trgt[sglx]];
								see_y = bgy[trgt[sglx]];
								see_z = bgz[trgt[sglx]];
								watch (); watch ();
								alfa += ashf[sglx];
								beta += bshf[sglx];
							}
						}
						break;
					case 'v':
						if (msg) {
							store_input_value ();
							subglxpar = 0;
							ccdpar = 0;
							glxpar = 2;
							sprintf (input, "%1.0f", vlct[sglx] / S);
							ic = strlen (input) - 1;
							point = 1;
							sign = 1;
							ics = 7;
						}
						else {
							glow = 1 - glow;
							if (glow) {
								fh = _open ("glowstar.pic", 0);
								if (fh==-1)
									glow = 0;
								else {
									initvesa (_640);
									if (video_mode[0]!=_640) {
									desp:	_640_480_16 ();
										_close (fh);
										glow = 0;
									}
									else {
										tavola_colori (gpal, 0, 256, 63, 63, 63);
										_read (fh, gws, 6594);
										_close (fh);
										fh = _creat ("vel_area.swp", 0);
										if (fh>-1) {
											_write (fh, vx, ms*g*8);
											_write (fh, vy, ms*g*8);
											i = _write (fh, vz, ms*g*8);
											if (i!=ms*g*8) goto desp;
											_close (fh);
											farfree (vx);
											farfree (vy);
											farfree (vz);
											farfree (adapted);
											desperation = (unsigned char far *) farmalloc (307232);
											if (!desperation) {
												vx = (double far *) farmalloc ((long)ms*(long)M_GAL*8);
												vy = (double far *) farmalloc ((long)ms*(long)M_GAL*8);
												vz = (double far *) farmalloc ((long)ms*(long)M_GAL*8);
												adapted = (unsigned char far *) farmalloc (38400);
												fh = _open ("vel_area.swp", 0);
												_read (fh, vx, ms*g*8);
												_read (fh, vy, ms*g*8);
												_read (fh, vz, ms*g*8);
												goto desp;
											}
											asm {	les ax, desperation
												mov seg_desperation, es
												inc seg_desperation }
											anchor = -1;
										}
										else
											goto desp;
									}
								}
							}
							else {
								farfree (desperation);
								vx = (double far *) farmalloc ((long)ms*(long)M_GAL*8);
								vy = (double far *) farmalloc ((long)ms*(long)M_GAL*8);
								vz = (double far *) farmalloc ((long)ms*(long)M_GAL*8);
								adapted = (unsigned char far *) farmalloc (38400);
								fh = _open ("vel_area.swp", 0);
								_read (fh, vx, ms*g*8);
								_read (fh, vy, ms*g*8);
								_read (fh, vz, ms*g*8);
								_close (fh);
								_640_480_16 ();
							}
						}
						break;
					case 'q':
						if (msg) {
							store_input_value ();
							subglxpar = 0;
							ccdpar = 0;
							glxpar = 5;
							sprintf (input, "%d", n[sglx]);
							ic = strlen (input) - 1;
							point = 0;
							sign = 0;
							ics = 4;
						}
						break;
					case '1':
						if (glxpar==4&&!subglxpar) {
							store_input_value ();
							sprintf (input, "%1.0f", (float)d[sglx]);
							ic = strlen (input) - 1;
							subglxpar = 1;
							ccdpar = 0;
							key = -1;
							ics = 6;
							sign = 0;
							point = 0;
						}
						break;
					case '2':
						if (glxpar==4&&!subglxpar) {
							store_input_value ();
							sprintf (input, "%1.2f", (float)(s[sglx] / 1E9));
							ic = strlen (input) - 1;
							subglxpar = 2;
							ccdpar = 0;
							key = -1;
							ics = 4;
							sign = 0;
							point = 1;
						}
						break;
					case 'c':
						if (!msg) {
							ccd_settings = 1 - ccd_settings;
							ccdpar = 0;
						}
						break;
					case 'k':
						if (msg) {
							store_input_value ();
							subglxpar = 0;
							glxpar = 1;
						}
						break;
					case 'i':
						if (msg) {
							store_input_value ();
							subglxpar = 0;
							glxpar = 6;
						}
						break;
					case 'e':
						if (ccd_settings) {
							store_input_value ();
							if (glxpar==2) glxpar = 0;
							if (subglxpar) subglxpar = 0;
							ccdpar = 1;
							sprintf (input, "%u", exposure);
							ic = strlen (input) - 1;
							point = 0;
							sign = 0;
							ics = 2;
						}
						break;
					case 'l':
						if (ccd_settings) {
							store_input_value ();
							if (glxpar==2) glxpar = 0;
							if (subglxpar) subglxpar = 0;
							ccdpar = 3;
							sprintf (input, "%d", pcrlim);
							ic = strlen (input) - 1;
							point = 0;
							sign = 0;
							ics = 2;
						}
						break;
					case 32:
						if (msg) {
							bgx[sglx] = cam_x;
							bgy[sglx] = cam_y;
							bgz[sglx] = cam_z;
							gx[sglx]  = cam_x;
							gy[sglx]  = cam_y;
							gz[sglx]  = cam_z;
							s_init (sglx, 0);
						}
						break;
					case 13:
						if (ccd_settings) {
							store_input_value ();
							ccd_image ();
							simwait ();
						}
						break;
					case '?':
						help = 1 - help;
						break;
					case 27:
						if (glxpar) {
							store_input_value ();
							key = -1;
							if (subglxpar)
								subglxpar = 0;
							else
								glxpar = 0;
							s_init (sglx, 0);
							break;
						}
						if (msg) {
							thrower ();
							g_init ();
							key = -1;
							msg = 0;
							break;
						}
						if (ccdpar) {
							store_input_value ();
							ccdpar = 0;
							key = -1;
							break;
						}
						if (lens) {
							lens = 0;
							key = -1;
							break;
						}
						if (ccd_settings) {
							ccd_settings = 0;
							key = -1;
							break;
						}
						if (gas) {
							gas = 0;
							key = -1;
							break;
						}
						if (help) {
							help = 0;
							key = -1;
							break;
						}
						if (key==27) exitflag = 1;
				}
				if ((key>='0'&&key<='9')||key=='.'||key=='-') {
					if (key=='-'&&(ic||input[0]!='0'||!sign)) goto denial;
					if (key=='.') {
						if (!point) goto denial;
						for (i=0; i<ic; i++)
							if (input[i]=='.')
								goto denial;
					}
					if (!ic&&input[ic]=='0'&&key!='.')
						ic = -1;
					else {
						j = 0;
						for (i=0; i<=ic; i++)
							if (input[i]!='.'&&input[i]!='-')
								j++;
						if (j==ics) goto denial;
					}
					ic++;
					input[ic] = key;
					input[ic+1] = 0;
					key = -1;
					denial:
				}
			}
			else {
				key = attendi_pressione_tasto();
				if (key==83) {
					if (!msg) goto lastone;
					if (g<2) goto lastone;
					if (sglx==anchor/ms) anchor = -1;
					if (sglx!=g-1) {
						j       = sglx;
						i       = sglx + 1;
					scroll:	t[j]    = t[i];
						s[j]    = s[i];
						d[j]    = d[i];
						n[j]    = n[i];
						gx[j]   = gx[i];
						gy[j]   = gy[i];
						gz[j]   = gz[i];
						gvx[j]  = gvx[i];
						gvy[j]  = gvy[i];
						gvz[j]  = gvz[i];
						bgx[j]  = bgx[i];
						bgy[j]  = bgy[i];
						bgz[j]  = bgz[i];
						vlct[j] = vlct[i];
						trgt[j] = trgt[i];
						ashf[j] = ashf[i];
						bshf[j] = bshf[i];
						_CX = n[i];
						_DI = ((ms * j) << 2) + 4;
						_SI = ((ms * i) << 2) + 4;
						asm {
							push ds
							mov es, seg_m
							mov ds, seg_m
							db 0xf3
							db 0x66
							db 0xa5
							pop ds
						}
						_CX = n[i] << 1;
						_DI = ((ms * j) << 3) + 4;
						_SI = ((ms * i) << 3) + 4;
						asm {
							push ds
							mov es, seg_x
							mov ds, seg_x
							db 0xf3
							db 0x66
							db 0xa5
							pop ds
						}
						_CX = n[i] << 1;
						_DI = ((ms * j) << 3) + 4;
						_SI = ((ms * i) << 3) + 4;
						asm {
							push ds
							mov es, seg_y
							mov ds, seg_y
							db 0xf3
							db 0x66
							db 0xa5
							pop ds
						}
						_CX = n[i] << 1;
						_DI = ((ms * j) << 3) + 4;
						_SI = ((ms * i) << 3) + 4;
						asm {
							push ds
							mov es, seg_z
							mov ds, seg_z
							db 0xf3
							db 0x66
							db 0xa5
							pop ds
						}
						_CX = n[i] << 1;
						_DI = ((ms * j) << 3) + 4;
						_SI = ((ms * i) << 3) + 4;
						asm {
							push ds
							mov es, seg_vx
							mov ds, seg_vx
							db 0xf3
							db 0x66
							db 0xa5
							pop ds
						}
						_CX = n[i] << 1;
						_DI = ((ms * j) << 3) + 4;
						_SI = ((ms * i) << 3) + 4;
						asm {
							push ds
							mov es, seg_vy
							mov ds, seg_vy
							db 0xf3
							db 0x66
							db 0xa5
							pop ds
						}
						_CX = n[i] << 1;
						_DI = ((ms * j) << 3) + 4;
						_SI = ((ms * i) << 3) + 4;
						asm {
							push ds
							mov es, seg_vz
							mov ds, seg_vz
							db 0xf3
							db 0x66
							db 0xa5
							pop ds
						}
						j++;
						i++;
						if (i<g) goto scroll;
					}
					g--;
					if (sglx) sglx--;
					msg = 6;
			lastone:}
				else {
					switch (key) {
						case 0x49: vrif_y -= 10*A; cam_y += vrif_y; break;
						case 0x51: vrif_y += 10*A; cam_y += vrif_y; break;
						case 0x8d:
							if (TT==0&&single_step<3*MIL) {
								single_step *= 2;
								S = 1000 * single_step;
								for (j=0; j<g; j++) vlct[j] *= 2;
								thrower ();
								g_init ();
							}
							break;
						case 0x91:
							if (TT==0&&single_step>0.2*MIL) {
								single_step /= 2;
								S = 1000 * single_step;
								for (j=0; j<g; j++) vlct[j] /= 2;
								thrower ();
								g_init ();
							}
							break;
						case 59:
							if (lens) {
								ob = 210;
								break;
							}
							if (gas) {
								if (file_selector("S|A|V|I|N|G| S|I|M|U|L|A|T|I|O|N|:", simfile, "SIM", SAVE)) {
									fh = _creat (file, 0);
									if (fh>-1) {
										_write (fh, &single_step, 8);
										_write (fh, &S, 	  8);
										_write (fh, &cam_x,      24);
										_write (fh, &alfa,       16);
										_write (fh, &delta_alfa, 16);
										_write (fh, &delta_pcam,  8);
										_write (fh, &ob,	 16);
										_write (fh, &anchor,	  2);
										_write (fh, &follow,      1);
										_write (fh, &g,		  2);
										_write (fh, n,		2*g);
										_write (fh, t,		2*g);
										_write (fh, X,		8*g);
										_write (fh, d,		8*g);
										_write (fh, s,		8*g);
										_write (fh, trgt,	2*g);
										_write (fh, ashf,	8*g);
										_write (fh, bshf,	8*g);
										_write (fh, vlct,	8*g);
										_write (fh, incl,	8*g);
										_write (fh, tilt,	8*g);
										_write (fh, acws,	  g);
										_write (fh, gx,		8*g);
										_write (fh, gy,		8*g);
										_write (fh, gz,		8*g);
										_write (fh, gvx,	8*g);
										_write (fh, gvy,	8*g);
										_write (fh, gvz,	8*g);
										_write (fh, bgx,	8*g);
										_write (fh, bgy,	8*g);
										_write (fh, bgz,	8*g);
										if (TT>0) {
											_write (fh, &TT, 8);
											for (j = 0; j < g; j++) {
												_write (fh, m+ms*j,  8*n[j]);
												_write (fh, x+ms*j,  8*n[j]);
												_write (fh, y+ms*j,  8*n[j]);
												_write (fh, z+ms*j,  8*n[j]);
												_write (fh, vx+ms*j, 8*n[j]);
												_write (fh, vy+ms*j, 8*n[j]);
												_write (fh, vz+ms*j, 8*n[j]);
											}
										}
										_close (fh);
										strcpy (simfile, file);
										pclear (adapted, 0);
										Testo (0, 470, "FILE STORED.", CENTRATO);
									}
									else
										Testo (0, 470, "UNABLE TO ACCESS DISK: DISK WRITE PROTECTED?", CENTRATO);
									QUADWORDS = 120;
									pcopy (adaptor+37600, adapted+37600);
									QUADWORDS = 9600;
									simwait ();
								}
							}
							break;
						case 60:
							if (lens) {
								ob = 300;
								break;
							}
							if (gas) {
								if (file_selector("L|O|A|D|I|N|G| S|I|M|U|L|A|T|I|O|N|:", simfile, "SIM", LOAD)) {
									fh = _open (file, 0);
									_read (fh, &single_step,  8);
									_read (fh, &S, 	  	  8);
									_read (fh, &cam_x,       24);
									_read (fh, &alfa,        16);
									_read (fh, &delta_alfa,  16);
									_read (fh, &delta_pcam,   8);
									_read (fh, &ob,	 	 16);
									_read (fh, &anchor,	  2);
									_read (fh, &follow,       1);
									_read (fh, &g,		  2);
									_read (fh, n,		2*g);
									_read (fh, t,		2*g);
									_read (fh, X,		8*g);
									_read (fh, d,		8*g);
									_read (fh, s,		8*g);
									_read (fh, trgt,	2*g);
									_read (fh, ashf,	8*g);
									_read (fh, bshf,	8*g);
									_read (fh, vlct,	8*g);
									_read (fh, incl,	8*g);
									_read (fh, tilt,	8*g);
									_read (fh, acws,	  g);
									_read (fh, gx,		8*g);
									_read (fh, gy,		8*g);
									_read (fh, gz,		8*g);
									_read (fh, gvx,		8*g);
									_read (fh, gvy,		8*g);
									_read (fh, gvz,		8*g);
									_read (fh, bgx,		8*g);
									_read (fh, bgy,		8*g);
									i = _read (fh, bgz,	8*g);
									if (!eof(fh)) {
										_read (fh, &TT, 8);
										for (j = 0; j < g; j++) {
											_read (fh, m+ms*j,  8*n[j]);
											_read (fh, x+ms*j,  8*n[j]);
											_read (fh, y+ms*j,  8*n[j]);
											_read (fh, z+ms*j,  8*n[j]);
											_read (fh, vx+ms*j, 8*n[j]);
											_read (fh, vy+ms*j, 8*n[j]);
											i = _read (fh, vz+ms*j, 8*n[j]);
										}
										if (i!=8*n[g-1]) {
											last_exit ();
											_close (fh);
											return;
										}
									}
									else {
										if (i!=8*g) {
											_close (fh);
											last_exit ();
											return;
										}
										else {
											TT = 0;
											thrower ();
											g_init ();
										}
									}
									_close (fh);
									gas = 0;
									msg = 0;
									lens = 0;
									help = 0;
									ccdpar = 0;
									glxpar = 0;
									subglxpar = 0;
									ccd_settings = 0;
									wait = 1;
									strcpy (simfile, file);
								}
							}
							break;
						case 61:
							if (lens) {
								ob = 810;
								break;
							}
							if (gas) {
								if (file_selector("D|E|L|E|T|I|N|G| S|I|M|U|L|A|T|I|O|N|:", simfile, "SIM", DELETE))
									remove (file);
							}
							break;
						case 62:
							if (gas) {
								if (file_selector("L|O|A|D|I|N|G S|I|M|U|L|A|T|E|D C|C|D| I|M|A|G|E|:", ccd_outfile, "BMP", LOAD)) {
									pclear (adapted, 0);
									fh = _open (file, 0);
									if (fh>-1) {
										_read (fh, tx, 54);
										if (memcmp(tx, bmpheader, 54)) {
											_close (fh);
											Testo (0, 470, "F|O|R|M|A|T| E|R|R|O|R|: THIS BITMAP IS NOT A STARWAYS' SIMULATED CCD PICTURE!", CENTRATO);
											goto wfmt;
										}
										grays (64);
										lseek (fh, 64, SEEK_CUR);
										vptr = 480; while (vptr) {
											_read (fh, adapted, 320);
											asm {	les si, adapted
												mov dx, vptr
												dec dx
												mov di, 640
												xor cx, cx
												xor bh, bh
												mov bl, 0x0F }
										criga:	asm {	mov al, es:[si]
												and al, bl
												not bl
												cmp bl, 0xF0
												jne stloc
												inc si
												jmp nxloc }
										stloc:  asm  	shr al, 4
										nxloc:	asm {	mov ah, 12
												int 0x10
												inc cx
												dec di
												jnz criga }
											vptr--;
										}
										_close(fh);
										Testo (0, 470, "IMAGE LOADED.", CENTRATO);
										strcpy (ccd_outfile, file);
									}
								wfmt:	QUADWORDS = 120;
									pcopy (adaptor+37600, adapted+37600);
									QUADWORDS = 9600;
									simwait ();
								}
							}
							break;
						case 63:
							if (gas) {
								if (file_selector("D|E|L|E|T|E| S|I|M|U|L|A|T|E|D| C|C|D I|M|A|G|E|:", ccd_outfile, "BMP", DELETE))
									remove (file);
							}
					}
				}
			}
			while (tasto_premuto())
				attendi_pressione_tasto();
		}
		else
			vrif_y = 0;

		/* Calcoli globali in risposta ai comandi impartiti. */

		if (fpcs) {
			k = M_PI / (200 * fpcs);
			w = (200 * A) / fpcs;
			q = 60 / fpcs;
		}
		else {
			k = M_PI / 5000;
			q = 25 / fpcs;
			w = 25 * A;
		}

		if (anchor==65535||!follow) {
			if (glxpar!=3&&glxpar!=6&&glxpar!=7) {
				delta_alfa -= (double) mdlty * k;
				delta_beta -= (double) mdltx * k;
				cam_x 	   += delta_xpos * q;
				cam_y 	   += delta_ypos * q;
				cam_z 	   += delta_zpos * q;
				alfa  	   += delta_alfa;
				beta  	   += delta_beta;
				if (mpul==1) {
					if (delta_pcam<0)
						delta_pcam = 0;
					else {
						if (delta_pcam||P!=1)
							delta_pcam += w;
					}
				}
				if (mpul==2) {
					if (delta_pcam>0)
						delta_pcam = 0;
					else {
						if (delta_pcam||P!=2)
							delta_pcam -= w;
					}
				}
				delta_xpos = - delta_pcam * sin(beta) * cos(alfa);
				delta_zpos =   delta_pcam * cos(beta) * cos(alfa);
				delta_ypos =   delta_pcam * sin(alfa);
				delta_alfa /= 1.05;
				delta_beta /= 1.05;
			}
			else {
				switch (glxpar) {
					case 3:
						cam_x       = bgx[sglx];
						cam_y       = bgy[sglx];
						cam_z       = bgz[sglx];
						delta_alfa  = 0;
						delta_beta  = 0;
						alfa       -= (double) mdlty * k;
						beta       -= (double) mdltx * k;
						delta_xpos  = 0;
						delta_ypos  = 0;
						delta_zpos  = 0;
						delta_pcam  = 0;
						break;
					case 6:
						alfa = 0;
						beta = M_PI_2;
						incl[sglx] -= (double) mdltx / DEG;
						s_init (sglx, 0);
						cam_x = bgx[sglx] + A/2 * d[sglx] * sin(beta) * cos(alfa);
						cam_z = bgz[sglx] - A/2 * d[sglx] * cos(beta) * cos(alfa);
						cam_y = bgy[sglx] - A/2 * d[sglx] * sin(alfa);
						break;
					case 7:
						beta = 0;
						alfa = M_PI_2;
						tilt[sglx] += (double) mdltx / DEG;
						s_init (sglx, 0);
						cam_x = bgx[sglx] + A * d[sglx] * sin(beta) * cos(alfa);
						cam_z = bgz[sglx] - A * d[sglx] * cos(beta) * cos(alfa);
						cam_y = bgy[sglx] - A * d[sglx] * sin(alfa);
				}
			}
		}

		if (anchor!=65535) {
			cam_x 	   =  x[anchor];
			cam_y 	   =  y[anchor];
			cam_z 	   =  z[anchor];
			delta_xpos = vx[anchor];
			delta_ypos = vy[anchor];
			delta_zpos = vz[anchor];
			if (follow==+1) {
				see_x = cam_x + (delta_xpos - gvx[anchor/ms]);
				see_y = cam_y + (delta_ypos - gvy[anchor/ms]);
				see_z = cam_z + (delta_zpos - gvz[anchor/ms]);
				watch ();
			}
			if (follow==-1) {
				see_x = cam_x - (delta_xpos - gvx[anchor/ms]);
				see_y = cam_y - (delta_ypos - gvy[anchor/ms]);
				see_z = cam_z - (delta_zpos - gvz[anchor/ms]);
				watch ();
			}
			if (follow==+2) {
				see_x = gx[anchor/ms];
				see_y = gy[anchor/ms];
				see_z = gz[anchor/ms];
				watch ();
			}
			if (wait) {
				delta_xpos = 0;
				delta_ypos = 0;
				delta_zpos = 0;
			}
		}

		rx      = cam_x - rif_x;
		ry      = cam_y - rif_y;
		rz      = cam_z - rif_z;
		cam_vel = sqrt (rx*rx + ry*ry + rz*rz);
		rif_x   = cam_x;
		rif_y   = cam_y;
		rif_z   = cam_z;

		/* Riempimento schermo. */

		if (!glow) pclear (adapted, 0);

		// Traccia il mirino d'inquadratura.

		if (!glow&&mpul==3) {
			asm mov es, seg_adapted;
			for (w=0; w<2*M_PI; w+=M_PI/50) {
				x2d = 75 * cos(w) + 320;
				shift = ((int)(75 * sin(w)) + 240) * 80 + 4;
				plot ();
			}
			for (w=0; w<2*M_PI; w+=M_PI/30) {
				x2d = 25 * cos(w) + 320;
				shift = ((int)(25 * sin(w)) + 240) * 80 + 4;
				plot ();
			}
			for (w=0; w<2*M_PI; w+=M_PI_4) {
				x2d = 27 * cos(w) + 320;
				shift = ((int)(27 * sin(w)) + 240) * 80 + 4;
				plot ();
			}
			for (w=0; w<2*M_PI; w+=M_PI_2) {
				x2d = 29 * cos(w) + 320;
				shift = ((int)(29 * sin(w)) + 240) * 80 + 4;
				plot ();
			}
		}

		if (glxpar==3) {
			Stick (300, 220, 340, 260);
			Stick (340, 220, 300, 260);
		}

		// Opzioni menó principale (sensibili al contesto).

		strcpy (escfunc, "E|S|C| @ END SESSION.");
		if (msg||ccd_settings||lens||gas) strcpy (escfunc, "E|S|C| @ PROGRESSIVELY CLOSES INTERFACE'S BUBBLES.");

		if (anchor!=65535) {
			i = 16;
			Testo (0, i, "D|ETACH FROM ANCHOR-STAR;", NORMALE);
		}
		else
			i = 8;

		Testo (0, i+8,  "M|ODIFY GALACTIC SCHEMES;", NORMALE);

		if (msg)
			i -= 8;
		else
			Testo (0, i+16, "R|EPEAT EXPERIMENT;", NORMALE);

		Testo (0, i+24, "O|BJECTIVE/LENS CHANGE;", NORMALE);

		if (!msg) {
			Testo (0, i+32, "A|RCHIVE ACCESS;", NORMALE);
			Testo (0, i+40, "C|CD SETUP;", NORMALE);
		}
		else
			i -= 16;

		if (grid)
			Testo (0, i+48, "G|RID O|N|/OFF;", NORMALE);
		else
			Testo (0, i+48, "G|RID ON/O|F|F|;", NORMALE);

		if (!wait) {
			Testo (0, i+56, "S|USPEND SIMULATION;", NORMALE);
			Testo (0, i+72, escfunc, NORMALE);
		}
		else {
			Testo (0, i+56, "P|ROCEED WITH SIMULATION;", NORMALE);
			Testo (0, i+64, "F|RAME BY FRAME;", NORMALE);
			if (TT>0)
				Testo (0, i+80, escfunc, NORMALE);
			else {
				if (clock()%5<3) {
					Testo (0, i+80, "C|T|R|L| + UP ARROW    R|A|I|S|E| TIME STEP.", NORMALE);
					Testo (0, i+88, "C|T|R|L| + DOWN ARROW  D|R|O|P| TIME STEP.", NORMALE);
				}
				Testo (0, i+104, escfunc, NORMALE);
			}
		}

		Testo (0, 450, "H|S|P|@S|T|A|R|W|A|Y|S|@1|9|9|7|-|2|0|0|3|E|N|", ALL_DX);

		pcs++;

		if (clock() - tci > 18) {
			if (!wait)
				sprintf (pcsline, "STEPS@PER@SECOND: %d", pcs);
			else
				sprintf (pcsline, "STEPS@PER@SECOND: 0 / P|A|U|S|E|D|");
			tci = clock ();
			if (pcs)
				fpcs = pcs;
			else
				fpcs = 1;
			pcs = 0;
		}

		Testo (0, 458, pcsline, ALL_DX);

		if (TT>0||clock()%5<3) {
			sprintf (tx, "YEARS@PER@STEP: %1.2f MILLION", single_step / MIL);
			Testo (0, 466, tx, ALL_DX);
		}

		rx = TT / (1000 * MIL);
		sprintf (tx, "TOTAL SIMULATION TIME: %1.1f BILLION YEARS", rx);
		Testo (0, 474, tx, ALL_DX);

		if (started) {
			if (rx>1.5) {
				started = 0;
				anchor = -1;
				follow = 0;
				help = 1;
			}
		}

		if (anchor==65535) {
			sprintf (tx, "FREE FLIGHT @ %04.0f KM/SEC", cam_vel / S);
			Testo (0, 474, tx, NORMALE);
		}
		else {
			Testo (0, 426, "N| @ NORMAL VIEW", NORMALE);
			Testo (0, 434, "X| @ CHASE FLIGHT DIRECTION", NORMALE);
			Testo (0, 442, "Y| @ CHASE OPPOSITE DIRECTION", NORMALE);
			Testo (0, 450, "Z| @ CHASE ORIGIN GALAXY", NORMALE);
			sprintf (tx, "CAMERA ATTACHED TO: STAR N.%d @ GALAXY N.%d", anchor%ms+1, anchor/ms+1);
			Testo (0, 466, tx, NORMALE);
			rx = vx[anchor] - gvx[anchor/ms];
			ry = vy[anchor] - gvy[anchor/ms];
			rz = vz[anchor] - gvz[anchor/ms];
			k = sqrt (rx*rx + ry*ry + rz*rz);
			sprintf (tx, "ACTUAL VELOCITY: %05.0f KM/SEC RELATIVE TO GALACTIC CORE", k / S);
			Testo (0, 474, tx, NORMALE);
			if (follow==+1) Testo (0, 0, "CHASING FLIGHT DIRECTION", ALL_DX);
			if (follow==-1) Testo (0, 0, "CHASING DIRECTION OPPOSITE TO FLIGHT", ALL_DX);
			if (follow==+2) {
				sprintf (tx, "CHASING ORIGIN GALAXY @ MASS=%1.2f BSM", (float)(s[anchor/ms] / 1E9));
				Testo (0, 0, tx, ALL_DX);
				rx = cam_x - gx[anchor/ms];
				ry = cam_y - gy[anchor/ms];
				rz = cam_z - gz[anchor/ms];
				f1 = rx*rx + ry*ry + rz*rz;
				if (f1<A) f1 = A;
				sprintf (tx, "DISTANCE FROM CENTER: %1.3f LY", (float)(sqrt(f1) / A));
				Testo (0, 16, tx, ALL_DX);
				w = (s[anchor/ms] * M * G * single_step * single_step) / f1;
				sprintf (tx, "GRAVITATIONAL PULL: %1.5f KM/SEC^", (float)(w / S));
				Testo (0, 8, tx, ALL_DX);
			}
		}

		cosbeta = cos(beta);
		sinbeta = sin(beta);
		cosalfa = cos(alfa);
		sinalfa = sin(alfa);

		cosbeta_x_ob = cos(beta) * ob;
		sinbeta_x_ob = sin(beta) * ob;
		cosalfa_x_ob = cos(alfa) * ob;
		sinalfa_x_ob = sin(alfa) * ob;

		if (glow) {
			if (hylit) {
				asm {   push ds
					mov ds, seg_desperation
					db 0x66; xor ax, ax
					mov si, 65532
					mov bx, ds
					xor dx, dx }
			clr4bs: asm {   db 0x66, 0x81, 0x24, 0xfe, 0xfe, 0xfe, 0xfe
					db 0x66; shr word ptr [si], 1
					sub si, 4
					jnc clr4bs
					add bx, 0x1000
					mov ds, bx
					inc dx
					cmp dx, 4
					ja  clrstp
					jne clr4bs
					mov si, 45052
					jmp clr4bs }
			clrstp: asm	pop ds
			}
			else {
				asm {   push ds
					mov ds, seg_desperation
					db 0x66; xor ax, ax
					mov si, 65532
					mov bx, ds
					xor dx, dx }
			_lr4bs: asm {   db 0x66; mov [si], ax
					sub si, 4
					jnc _lr4bs
					add bx, 0x1000
					mov ds, bx
					inc dx
					cmp dx, 4
					ja  _lrstp
					jne _lr4bs
					mov si, 45052
					jmp _lr4bs }
			_lrstp: asm	pop ds
			}
			glowin_render ();
			pselect (0);
			asm {   push ds
				les di, dword ptr adaptor
				mov ds, seg_desperation
				mov cx, 16384
				mov bx, ds
				xor si, si
				xor dx, dx }
		copypg: asm {   db 0xf3 // macro: rep movsd
				db 0x66
				db 0xa5
				mov cx, 16384
				add bx, 0x1000
				mov ds, bx
				inc dx
				push ax
				push bx
				mov ax, 0x4f05
				mov bx, 1
				int 0x10
				pop bx
				pop ax
				cmp dx, 4
				ja  copyst
				jne copypg
				mov cx, 11263
				jmp copypg }
		copyst: asm	pop ds
		}
		else
			render ();

		if (grid&&!glow) {

			q = 1E5*A;
			k = q / 2;
			w = -14*q;
			ww = 14*q;

			airbrush = 1;

			stick3d (cam_x, 0, w, cam_x, 0, ww);
			stick3d (w, 0, cam_z, ww, 0, cam_z);
			stick3d (w, cam_y, 0, ww, cam_y, 0);
			stick3d (0, cam_y, w, 0, cam_y, ww);

			for (i=-14, f1=w; f1<=0; f1+=q, i++) {

				airbrush = abs(i) + 2;

				stick3d (f1, 0, w, f1, 0, ww);
				stick3d (w, 0, f1, ww, 0, f1);
				stick3d (f1, w, 0, f1, ww, 0);
				stick3d (0, w, f1, 0, ww, f1);

				stick3d (-f1, 0, w, -f1, 0, ww);
				stick3d (w, 0, -f1, ww, 0, -f1);
				stick3d (-f1, w, 0, -f1, ww, 0);
				stick3d (0, w, -f1, 0, ww, -f1);

				if (airbrush>1) {

					stick3d (0, f1, w, 0, f1, ww);
					stick3d (w, f1, 0, ww, f1, 0);

					stick3d (0, -f1, w, 0, -f1, ww);
					stick3d (w, -f1, 0, ww, -f1, 0);

				}

			}

		}

		/* Elementi in sovraimpressione. */

		if (msg) {
			Testo (350, 285, "S|P|A|C|E| @ MOVES SELECTED GALAXY.", NORMALE);
			i = 315;
			Testo (350, 300, "D|I|S|T|A|N|C|E| F|R|O|M| C|A|M|E|R|A|:", NORMALE);
			for (j=0; j<g; j++) {
				f1 = cam_x - bgx[j];
				f2 = cam_y - bgy[j];
				f3 = cam_z - bgz[j];
				f1 = sqrt (f1*f1+f2*f2+f3*f3);
				sprintf (tx, "FROM GALAXY N. %d| = %1.0f L|.Y|.", j+1, f1/A);
				Testo (365, i, tx, NORMALE);
				if (j==trgt[sglx]&&j!=sglx) {
					i += 8;
					mem ();
					cam_x = bgx[sglx];
					cam_y = bgy[sglx];
					cam_z = bgz[sglx];
					f1 = cam_x - bgx[j];
					f2 = cam_y - bgy[j];
					f3 = cam_z - bgz[j];
					f4 = sqrt (f1*f1+f2*f2+f3*f3);
					f5 = f4;
					see_x = bgx[j];
					see_y = bgy[j];
					see_z = bgz[j];
					watch ();
					watch ();
					alfa += ashf[sglx];
					beta += bshf[sglx];
					forward (f4);
					f1 = cam_x - bgx[j];
					f2 = cam_y - bgy[j];
					f3 = cam_z - bgz[j];
					f1 = sqrt (f1*f1+f2*f2+f3*f3);
					if (f1<f5) f5 = f1;
					forward (-2*f4);
					f1 = cam_x - bgx[j];
					f2 = cam_y - bgy[j];
					f3 = cam_z - bgz[j];
					f1 = sqrt (f1*f1+f2*f2+f3*f3);
					if (f1<f5) f5 = f1;
					rcm ();
					sprintf (tx, "TARGET NEAREST FLY-BY %1.0f L.Y.", f5/A);
					Testo (380, i, tx, NORMALE);
					i += 8;
					if (vlct[sglx]>S) {
						f1 = (f4/1000) / (vlct[sglx]/S); // Secondi all'avvicinamento.
						f1 /= 1000*MIL; // Miliardi di anni all'avvicinamento.
						sprintf (tx, "FLY-BY TIME EST. = %1.3f BILLION YRS", f1);
						Testo (380, i, tx, NORMALE);
						i += 8;
					}
					if (f4>100*A) {
						f1 = sqrt ((M * G * s[j]) / f4) / 1000;
						sprintf (tx, "ORBITAL VELOCITY = %1.3f KM/S", f1);
						Testo (380, i, tx, NORMALE);
						i += 8;
					}
				}
				i += 10;
			}
			bolla (130, 208, 125);
			Testo (-190, 120, "G|A|L|A|C|T|I|C| S|C|H|E|M|E|S|", CENTRATO);
			Testo (-190, 136, "TO SELECT GALAXY TO BE MODIFIED,", CENTRATO);
			Testo (-190, 144, "PRESS BOTH MOUSE BUTTONS, MOVE", CENTRATO);
			Testo (-190, 152, "THE SIGHT UNTIL THE SAID GALAXY'S", CENTRATO);
			Testo (-190, 160, "CORE LAYS IN THE SIGHT'S AREA,", CENTRATO);
			Testo (-190, 168, "THEN RELEASE THE MOUSE BUTTONS.", CENTRATO);
			sprintf (tx, "ACTUALLY SELECTED GALAXY IS N. %d:", sglx + 1);
			Testo (-190, 208, tx, CENTRATO);
			if (msg==1) sprintf (tx, "THE GALAXY WAS PREVIOUSLY");
			if (msg==2) sprintf (tx, "IT'S BEEN SELECTED WITH THE SIGHT.");
			if (msg==3) sprintf (tx, "THE AREA HIGHLIGHTED BY THE SIGHT");
			if (msg==4) sprintf (tx, "IT HAS BEEN JUST ADDED.");
			if (msg==5) sprintf (tx, "T|H|I|S| V|E|R|S|I|O|N| O|F| S|T|A|R|W|A|Y|S|");
			if (msg==6) sprintf (tx, "THE SELECTION HAS IN REALITY CHANGED");
			Testo (-190, 216, tx, CENTRATO);
			if (msg!=2&&msg!=4) {
				if (msg==1) sprintf (tx, "SELECTED (BY DEFAULT, THE FIRST ONE).");
				if (msg==3) sprintf (tx, "DIDN'T INCLUDE GALACTIC CORES.");
				if (msg==5) sprintf (tx, "D|O|E|S| N|O|T| M|A|N|A|G|E| M|O|R|E| T|H|A|N| 4| G|A|L|A|X|I|E|S|.|");
				if (msg==6) sprintf (tx, "BECAUSE SELECTED GALAXY WAS REMOVED.");
				Testo (-190, 224, tx, CENTRATO);
			}
			Testo (-190, 240, "PRESS +| TO A|D|D| A GALAXY,", CENTRATO);
			Testo (-190, 248, "D|E|L| TO D|E|L|E|T|E| SELECTED GALAXY.", CENTRATO);
			Testo (-190, 264, "ONCE YOU'RE DONE WITH CHANGES", CENTRATO);
			Testo (-190, 272, "PRESS E|S|C| TO CLOSE THIS", CENTRATO);
			Testo (-190, 280, "BUBBLE.", CENTRATO);
			bolla (280, 85, 80);
			Testo (-40, 20, "A|C|T|U|A|L|", CENTRATO);
			Testo (-40, 28, "S|E|L|E|C|T|I|O|N|:", CENTRATO);
			sprintf (tx, "GALAXY N. %d|", sglx + 1);
			Testo (-40, 36, tx, CENTRATO);
			Testo (-40, 46, "C|U|R|R|E|N|T| P|A|R|A|M|E|T|E|R|S|", CENTRATO);
			sprintf (tx, "K|LASS: %s", gtyp[t[sglx]]);
			Testo (212, 62, tx, NORMALE);
			sprintf (tx, "V|ELOCITY: %1.0f KM/S", vlct[sglx] / S);
			Testo (212, 70, tx, NORMALE);
			if (trgt[sglx]==sglx||trgt[sglx]>=g)
				sprintf (tx, "T|ARGET: N|O|N|E|");
			else
				sprintf (tx, "T|ARGET: GALAXY N.%d|", trgt[sglx] + 1);
			Testo (212, 78, tx, NORMALE);
			sprintf (tx, "I|NCLINATION: %1.1f\\", incl[sglx]*DEG);
			Testo (212, 86, tx, NORMALE);
			sprintf (tx, "A|ZIMUT: %1.1f\\", tilt[sglx]*DEG);
			Testo (212, 94, tx, NORMALE);
			sprintf (tx, "Q|UANTITY OF STARS: %d", n[sglx]);
			Testo (212, 102, tx, NORMALE);
			Testo (212, 114, "D|IMENSIONS:", NORMALE);
			sprintf (tx, "DIAMETER %1.0f LY", (float)d[sglx]);
			Testo (222, 122, tx, NORMALE);
			sprintf (tx, "MASS %1.2f BSM", (float)(s[sglx] / 1E9));
			Testo (222, 130, tx, NORMALE);
			if (glxpar) {
				switch (glxpar) {
					case 1:
						if (key==72&&t[sglx]>0) {
							t[sglx]--;
							s_init (sglx, 0);
						}
						if (key==80&&t[sglx]<3) {
							t[sglx]++;
							s_init (sglx, 0);
						}
						bolla (410, 80, 70);
						for (j = 50, i = 0; i < 4; i++, j += 10) {
							Testo (90, j, gtyp[i], CENTRATO);
							if (i==t[sglx]) Testo (91, j, gtyp[i], CENTRATO);
						}
						Testo (90, 105, "R|OTATION:", CENTRATO);
						if (acws[sglx])
							Testo (90, 118, "A|N|T|I|CLOCKWISE.", CENTRATO);
						else
							Testo (90, 118, "CLOCKWISE.", CENTRATO);
						break;
					case 2:
						bolla (410, 88, 70);
						Testo (90, 64, "@ V|E|L|O|C|I|T|Y| @", CENTRATO);
						sprintf (tx, "KM/SEC: %s|", input);
						Testo (90, 96, tx, CENTRATO);
						break;
					case 3:
						if (key==72&&trgt[sglx]>0) {
							trgt[sglx]--;
							if (trgt[sglx]!=sglx) {
								cam_x = bgx[sglx];
								cam_y = bgy[sglx];
								cam_z = bgz[sglx];
								see_x = bgx[trgt[sglx]];
								see_y = bgy[trgt[sglx]];
								see_z = bgz[trgt[sglx]];
								watch ();
								watch ();
							}
						}
						if (key==80&&trgt[sglx]<g-1) {
							trgt[sglx]++;
							if (trgt[sglx]!=sglx) {
								cam_x = bgx[sglx];
								cam_y = bgy[sglx];
								cam_z = bgz[sglx];
								see_x = bgx[trgt[sglx]];
								see_y = bgy[trgt[sglx]];
								see_z = bgz[trgt[sglx]];
								watch ();
								watch ();
							}
						}
						bolla (410, 96, 70);
						for (j = 64, i = 0; i < 4; i++, j += 12) {
							if (i==sglx)
								sprintf (tx, "NONE");
							else {
								if (i<g)
									sprintf (tx, "GALAXY N. %d", i + 1);
								else
									sprintf (tx, "(NOT APPLICABLE)");
							}
							Testo (90, j, tx, CENTRATO);
							if (i==trgt[sglx]) Testo (91, j, tx, CENTRATO);
						}
						if (trgt[sglx]==sglx) {
							ashf[sglx] = 0;
							bshf[sglx] = 0;
						}
						else {
							_alfa = alfa;
							_beta = beta;
							see_x = bgx[trgt[sglx]];
							see_y = bgy[trgt[sglx]];
							see_z = bgz[trgt[sglx]];
							watch ();
							watch ();
							ashf[sglx] = _alfa - alfa;
							bshf[sglx] = _beta - beta;
							alfa = _alfa;
							beta = _beta;
						}
						k = fabs (bshf[sglx] * DEG);
						if (k>360) k -= 360;
						if (k>180) k -= 360;
						sprintf (tx, "DELTA/R.A. = %1.1f\\", (float)fabs(k));
						Testo (90, j + 4, tx, CENTRATO);
						k = fabs (ashf[sglx] * DEG);
						if (k>360) k -= 360;
						if (k>180) k -= 360;
						sprintf (tx, "DELTA/DECL = %1.1f\\", (float)fabs(k));
						Testo (90, j + 12, tx, CENTRATO);
						break;
					case 4:
						bolla (410, 120, 70);
						Testo (90, 90, "D|I|M|E|N|S|I|O|N|S|", CENTRATO);
						Testo (350, 120, "1| @ DIAMETER", NORMALE);
						Testo (350, 130, "2| @ ATTRACTOR'S MASS", NORMALE);
						if (subglxpar==1) {
							bolla (500, 120, 40);
							Testo (180, 110, "D|I|A|M|E|T|E|R", CENTRATO);
							sprintf (tx, "%s| L.Y.", input);
							Testo (180, 130, tx, CENTRATO);
						}
						if (subglxpar==2) {
							bolla (500, 130, 40);
							Testo (180, 110, "C|O|R|E|", CENTRATO);
							Testo (180, 118, "M|A|S|S|", CENTRATO);
							Testo (180, 130, "(BSM)", CENTRATO);
							sprintf (tx, "%s|", input);
							Testo (180, 146, tx, CENTRATO);
						}
						break;
					case 5:
						bolla (410, 104, 70);
						Testo (90, 80, "N|R|.| O|F| S|T|A|R|S|", CENTRATO);
						Testo (90, 90, "(MAX. 1500)", CENTRATO);
						sprintf (tx, "%s|", input);
						Testo (90, 128, tx, CENTRATO);
						break;
					case 6:
						bolla (410, 112, 70);
						Testo (90, 88, "I|N|C|L|I|N|A|T|I|O|N|", CENTRATO);
						sprintf (tx, "%1.1f\\", incl[sglx]*DEG);
						Testo (90, 128, tx, CENTRATO);
						break;
					case 7:
						bolla (410, 120, 70);
						Testo (90, 96, "A|Z|I|M|U|T|", CENTRATO);
						sprintf (tx, "%1.1f\\", tilt[sglx]*DEG);
						Testo (90, 136, tx, CENTRATO);
						break;
				}
			}
		}

		if (ccd_settings) {
			bolla (150, 350, 85);
			Testo (-170, 302, "C|C|D S|E|T|U|P|", CENTRATO);
			sprintf (tx, "E|XPOSURE: %d MIN.", exposure);
			Testo (92, 326, tx, NORMALE);
			sprintf (tx, "F|OCUS DEVIATION: %d", (int) mag_factor);
			Testo (92, 334, tx, NORMALE);
			sprintf (tx, "L|IMIT DIFFRACTION: %d%", pcrlim);
			Testo (92, 342, tx, NORMALE);
			Testo (92, 350, "O|UTPUT TO FILE...", NORMALE);
			Testo (-170, 370, "TO OBTAIN THE IMAGE,", CENTRATO);
			Testo (-170, 378, "PRESS E|N|T|E|R|.", CENTRATO);
			Testo (-170, 390, "TO REMOVE THIS", CENTRATO);
			Testo (-170, 398, "BUBBLE, PRESS C|.", CENTRATO);
			switch (ccdpar) {
				case 1:
					bolla (270, 334, 60);
					Testo (-50, 316, "E|X|P|O|S|U|R|E|", CENTRATO);
					sprintf (tx, "MIN: %s|", input);
					Testo (-50, 348, tx, CENTRATO);
					break;
				case 2:
					bolla (270, 340, 60);
					Testo (-50, 316, "D|E|V|I|A|T|I|O|N|", CENTRATO);
					Testo (-50, 324, "F|R|O|M|", CENTRATO);
					Testo (-50, 332, "P|E|R|F|E|C|T| F|O|C|U|S|", CENTRATO);
					sprintf (tx, "%s|", input);
					Testo (-50, 356, tx, CENTRATO);
					break;
				case 3:
					bolla (270, 358, 60);
					Testo (-50, 332, "P|E|R|C|E|N|T|A|G|E|", CENTRATO);
					Testo (-50, 340, "O|F| S|T|A|R|S|", CENTRATO);
					Testo (-50, 348, "C|A|U|S|I|N|G|", CENTRATO);
					Testo (-50, 356, "D|I|F|F|R|A|C|T|I|O|N|:", CENTRATO);
					sprintf (tx, "%s|", input);
					Testo (-50, 380, tx, CENTRATO);
					break;
			}
		}

		if (lens) {
			bolla (540, 250, 90);
			Testo (220, 214, "C|H|O|O|S|E|", CENTRATO);
			Testo (220, 222, "C|A|M|E|R|A|", CENTRATO);
			Testo (220, 230, "O|B|J|E|C|T|I|V|E|/|L|E|N|S|", CENTRATO);
			if (ob==210)
				Testo (460, 258, "F|1| @  35 MM - W|I|D|E|-|A|N|G|L|E|", NORMALE);
			else
				Testo (460, 258, "F|1| @  35 MM - WIDE-ANGLE", NORMALE);
			if (ob==300)
				Testo (460, 266, "F|2| @  50 MM - N|O|R|M|A|L|", NORMALE);
			else
				Testo (460, 266, "F|2| @  50 MM - NORMAL", NORMALE);
			if (ob==810)
				Testo (460, 274, "F|3| @ 135 MM - T|E|L|E|", NORMALE);
			else
				Testo (460, 274, "F|3| @ 135 MM - TELE", NORMALE);
		}

		if (gas) {
			bolla (400, 300, 100);
			Testo (80, 260, "- S|I|M|U|L|A|T|I|O|N|S|'|A|R|C|H|I|V|E| -", CENTRATO);
			Testo (335, 280, "F|1| SAVE   SIMULATION.", NORMALE);
			Testo (335, 288, "F|2| LOAD   SIMULATION.", NORMALE);
			Testo (335, 296, "F|3| DELETE SIMULATION.", NORMALE);
			Testo (335, 312, "F|4| LOAD   PICTURE.", NORMALE);
			Testo (335, 320, "F|5| DELETE PICTURE.", NORMALE);
			Testo (80, 340, "A|: CLOSE THIS BUBBLE.", CENTRATO);
		}

		if (exitflag) {
			bolla (320, 240, 200);
			Testo (0, 210, "E|N|D| S|E|S|S|I|O|N|!", CENTRATO);
			Testo (0, 240, "PRESS:", CENTRATO);
			Testo (272, 267, "Y| TO CONFIRM.", NORMALE);
			Testo (272, 280, "N| TO CANCEL.", NORMALE);
		}

		if (!help)
			Testo (0, 0, "?| @ O|N|L|I|N|E| H|E|L|P| (ACTUALLY DISABLED).", NORMALE);
		else {
			Testo (0, 0, "?| @ DISABLE ONLINE HELP BUBBLE.", NORMALE);
			if (msg) {
				bolla (500, 233, 120);
				Testo (180, 148, "@ H|E|L|P| @", CENTRATO);
				if (!glxpar) {
					Testo (180, 164, "IN THE BUBBLE TO THE LEFT", CENTRATO);
					Testo (180, 172, "IS INDICATED HOW TO SELECT,", CENTRATO);
					Testo (180, 180, "ADD, AND REMOVE GALAXIES.", CENTRATO);
					Testo (180, 196, "TO MODIFY PARAMETERS OF", CENTRATO);
					Testo (180, 204, "SELECTED GALAXY, SHOWEN IN", CENTRATO);
					Testo (180, 212, "THE SMALLER BUBBLE, PRESS", CENTRATO);
					Testo (180, 220, "ONE OF THE LETTERS HIGHLIGHTED", CENTRATO);
					Testo (180, 228, "IN THE SMALL BUBBLE.", CENTRATO);
					Testo (180, 244, "ADDED GALAXIES WILL BE PLACED", CENTRATO);
					Testo (180, 252, "IN THE POINT WHERE THE CAMERA", CENTRATO);
					Testo (180, 260, "IS ACTUALLY POSITIONED.", CENTRATO);
					Testo (180, 276, "IF THERE IS ONLY ONE GALAXY,", CENTRATO);
					Testo (180, 284, "IT'S IMPOSSIBLE TO REMOVE IT.", CENTRATO);
				}
				else {
					if (glxpar!=3) {
						switch (glxpar) {
							case 1:
								Testo (180, 164, "SELECT AN", CENTRATO);
								Testo (180, 172, "A|R|C|H|E|T|Y|P|E| OF GALAXY", CENTRATO);
								Testo (180, 180, "THROUGHT THE C|U|R|S|O|R| A|R|R|O|W|S|", CENTRATO);
								Testo (180, 188, "(UP ARROW, DOWN ARROW).", CENTRATO);
								Testo (180, 212, "KEY R| ALLOWS TO SPECIFY", CENTRATO);
								Testo (180, 220, "THE VERSE IN WHICH STARS ROTATE", CENTRATO);
								Testo (180, 228, "IN RESPECT TO A VIEW POINT PLACED", CENTRATO);
								Testo (180, 236, "RIGHT ABOVE THE GALAXY'S CORE.", CENTRATO);
								break;
							case 2:
								Testo (180, 164, "THIS PARAMETER", CENTRATO);
								Testo (180, 172, "REPRESENTS THE V|E|L|O|C|I|T|Y|", CENTRATO);
								Testo (180, 180, "AT WHICH THE GALAXY WILL MOVE", CENTRATO);
								Testo (180, 188, "TOWARD THE SPECIFIED TARGET.", CENTRATO);
								Testo (180, 204, "BELOW TO THE RIGHT, ON THE SCREEN,", CENTRATO);
								Testo (180, 212, "IS AN ESTIMATE OF THE TIME BEFORE", CENTRATO);
								Testo (180, 220, "THE GALAXY REACHES THE MINIMUM", CENTRATO);
								Testo (180, 228, "DISTANCE FROM ITS TARGET,", CENTRATO);
								Testo (180, 236, "AND THE SPEED NECESSARY", CENTRATO);
								Testo (180, 244, "FOR A CIRCULAR ORBIT.", CENTRATO);
								break;
							default:
								Testo (180, 164, "WRITE THE N|E|W| V|A|L|U|E", CENTRATO);
								Testo (180, 172, "FOR THE PARAMETER.", CENTRATO);
						}
						if (glxpar>2) {
							Testo (180, 204, "TO MODIFY THE VALUE,", CENTRATO);
							Testo (180, 212, "USE KEYS 0 TO 9, USE", CENTRATO);
							Testo (180, 220, "THE D|O|T| (.) TO MARK DECIMALS,", CENTRATO);
							Testo (180, 228, "AND H|Y|P|E|N| (-) TO MARK NEGATIVES.", CENTRATO);
							Testo (180, 236, "B|A|C|K|S|P|A|C|E| DELETES", CENTRATO);
							Testo (180, 244, "DIGITS STARTING FROM R|I|G|H|T|.", CENTRATO);
						}
						Testo (180, 260, "BEFORE CHANGING ANOTHER VALUE,", CENTRATO);
						Testo (180, 268, "MAKE SURE THE ACTUAL PARAMETER'S", CENTRATO);
						Testo (180, 276, "VALUE IS SET TO THE DESIRED VALUE.", CENTRATO);
						Testo (180, 300, "TO CLOSE THE BUBBLE FOR", CENTRATO);
						Testo (180, 308, "A CERTAIN PARAMETER", CENTRATO);
						Testo (180, 316, "PRESS E|S|C|.", CENTRATO);
					}
					else {
						Testo (180, 164, "THE CAMERA HAS BEEN", CENTRATO);
						Testo (180, 172, "PLACED IN THE MIDDLE OF THE", CENTRATO);
						Testo (180, 180, "CURRENTLY SELECTED GALAXY,", CENTRATO);
						Testo (180, 188, "SO THAT IT'S EASIER TO POINT", CENTRATO);
						Testo (180, 196, "THE GALAXY TO THE DIRECTION", CENTRATO);
						Testo (180, 204, "WHICH THE CAMERA ITSELF POINTS.", CENTRATO);
						Testo (180, 220, "ONCE THE TARGET HAS BEEN CHOSEN,", CENTRATO);
						Testo (180, 228, "USING CURSOR UP/DOWN KEYS, IT", CENTRATO);
						Testo (180, 236, "WILL BE POSSIBLE TO SPECIFY,", CENTRATO);
						Testo (180, 244, "QUITE PRECISELY, THE POINT WHERE", CENTRATO);
						Testo (180, 252, "THE FLY-BY WILL OCCUR, BY MOVING", CENTRATO);
						Testo (180, 260, "THE MOUSE AND WITH THE HELP OF THE", CENTRATO);
						Testo (180, 268, "CROSSHAIR.", CENTRATO);
						Testo (180, 284, "T|O| M|A|K|E| A|N| O|R|B|I|T|, SET", CENTRATO);
						Testo (180, 292, "THE D|E|L|T|A|/R|.A|. TO 90\\", CENTRATO);
						Testo (180, 300, "(HORIZONTAL PLANE ORBIT)", CENTRATO);
						Testo (180, 308, "OR D|E|L|T|A|/D|E|C|L|. TO 90\\", CENTRATO);
						Testo (180, 316, "(VERTICAL PLANE ORBIT).", CENTRATO);
					}
				}
			}
			else {
				bolla (500, 143, 120);
				Testo (180,  48, "@ H|E|L|P| @", CENTRATO);
				if (anchor==65535) {
					if (mpul!=3) {
						Testo (180,  64, "WELCOME TO S|T|A|R|W|A|Y|S|.", CENTRATO);
						Testo (180,  72, "THE MAIN OPTIONS ARE SHOWEN", CENTRATO);
						Testo (180,  80, "IN THE UPPER LEFT CORNER OF THE", CENTRATO);
						Testo (180,  88, "SCREEN. TO ACCESS ONE OF THEM,", CENTRATO);
						Testo (180,  96, "PRESS THE LETTER KEY CORRESPONDING", CENTRATO);
						Testo (180, 104, "TO THE H|I|G|H|L|I|G|H|T|E|D| LETTER OF", CENTRATO);
						Testo (180, 112, "THE OPTION'S NAME.", CENTRATO);
						Testo (180, 128, "TO MOVE THE CAMERA FORWARD AND BACK,", CENTRATO);
						Testo (180, 136, "PRESS THE L|E|F|T| MOUSE BUTTON OR", CENTRATO);
						Testo (180, 144, "THE R|I|G|H|T| MOUSE BUTTON; TOGGLING", CENTRATO);
						Testo (180, 152, "BUTTONS WILL STOP THE CAMERA.", CENTRATO);
						Testo (180, 160, "TO CHANGE DIRECTION, MOVE THE", CENTRATO);
						Testo (180, 168, "MOUSE; H|O|L|D|I|N|G| D|O|W|N|", CENTRATO);
						Testo (180, 176, "BOTH MOUSE BUTTONS, A SPECIAL", CENTRATO);
						Testo (180, 184, "SIGHT WILL APPEAR.", CENTRATO);
						Testo (180, 208, "I|F| T|H|I|S| B|U|B|B|L|E| COVERS", CENTRATO);
						Testo (180, 216, "SOME OF THE DATA, USE ?| TO", CENTRATO);
						Testo (180, 224, "TOGGLE IT ON/OFF.", CENTRATO);
					}
					else {
						Testo (180,  64, "WHEN BUTTONS WILL BE", CENTRATO);
						Testo (180,  72, "RELEASED, THE NEAREST", CENTRATO);
						Testo (180,  80, "STAR VISIBLE IN THE", CENTRATO);
						Testo (180,  88, "SMALLEST SIGHT'S CIRCLE", CENTRATO);
						Testo (180,  96, "WILL BECOME THE NEW", CENTRATO);
						Testo (180, 104, "VIEW POINT.", CENTRATO);
						Testo (180, 120, "THE CAMERA WILL THEN MOVE", CENTRATO);
						Testo (180, 128, "ALONG WITH THAT STAR.", CENTRATO);
					}
				}
				else {
					Testo (180,  64, "THE CAMERA IS NOW", CENTRATO);
					Testo (180,  72, "ATTACHED TO THE STAR WHOSE", CENTRATO);
					Testo (180,  80, "INFORMATIONS ARE VISIBLE IN", CENTRATO);
					Testo (180,  88, "THE LOWER LEFT CORNER.", CENTRATO);
					Testo (180, 104, "TO GET BACK TO FREE FLIGHT", CENTRATO);
					Testo (180, 112, "MODE, PRESS THE D| KEY.", CENTRATO);
					Testo (180, 128, "TO FOLLOW THE DIRECTION", CENTRATO);
					Testo (180, 136, "IN WHICH THE STAR IS MOVING,", CENTRATO);
					Testo (180, 144, "PRESS THE X| KEY.", CENTRATO);
					Testo (180, 160, "PRESSING Y| THE CAMERA WILL", CENTRATO);
					Testo (180, 168, "CHASE THE DIRECTION FROM", CENTRATO);
					Testo (180, 176, "WHICH THE STAR COMES FROM, WHILE", CENTRATO);
					Testo (180, 184, "THE Z| KEY WILL CHASE THE CORE", CENTRATO);
					Testo (180, 192, "OF THE GALAXY THE STAR IS", CENTRATO);
					Testo (180, 200, "PART OF.", CENTRATO);
					Testo (180, 216, "TO DISABLE X|, Y| OR Z|,", CENTRATO);
					Testo (180, 224, "PRESS THE N| KEY.", CENTRATO);
				}
			}
		}

		/* Controllo collisioni tra nuclei galattici. */

		if (g>1) {
			spx = 250;
			for (pan = 0; pan < g; pan++) { // Controllo velocitÖ luminale.
				w  = gvx[pan]*gvx[pan];
				w += gvy[pan]*gvy[pan];
				w += gvz[pan]*gvz[pan];
				w  = sqrtl (w) / S;
				if (w>=C||pan==errorflag) {
					bolla (320 + spx, 380, 55);
// Linea nr. 4444, carattere nr.44 ->      x   (carino, no?)
					Testo (spx, 361, "GALAXY NUMBER", CENTRATO);
					sprintf (tx, "%d|", pan+1);
					Testo (spx, 373, tx, CENTRATO);
					Testo (spx, 385, "HAS CROSSED C|!", CENTRATO);
					Testo (spx, 397, "(PRECISION LOSS)", CENTRATO);
					errorflag = pan;
					spx -= 115; if (spx<-250) goto underflow;
					wait = 1;
				}
			}
			for (pan = 0; pan < g; pan ++) { // pan = galassia bersaglio.
				for (aux = 0; aux < g; aux ++) { // aux = galassia proiettile.
					if (pan!=aux) {
						ww = gx[pan] - gx[aux]; ww *= ww;
						q  = gy[pan] - gy[aux]; ww += q * q;
						q  = gz[pan] - gz[aux]; ww += q * q;
						if (!ww) {
							if (aux>pan) {
								bolla (320 + spx, 380, 55);
								Testo (spx, 361, "GALAXIES N.", CENTRATO);
								sprintf (tx, "%d| AND %d|", pan+1, aux+1);
								Testo (spx, 373, tx, CENTRATO);
								Testo (spx, 385, "HAVE C|O|L|L|I|D|E|D|", CENTRATO);
								Testo (spx, 397, "(MELTED TOGETHER)", CENTRATO);
								spx -= 115; if (spx<-250) goto underflow;
							}
						}
						else {
							k = fabs(gvx[pan]) + fabs(gvx[aux]); k *= k;
							q = fabs(gvy[pan]) + fabs(gvy[aux]); k += q * q;
							q = fabs(gvz[pan]) + fabs(gvz[aux]); k += q * q;
							if (ww<=k) {
								w = (gx[pan] + gx[aux]) / 2;
								gx[pan] = w; gx[aux] = w;
								w = (gy[pan] + gy[aux]) / 2;
								gy[pan] = w; gy[aux] = w;
								w = (gz[pan] + gz[aux]) / 2;
								gz[pan] = w; gz[aux] = w;
								q = s[aux] / s[pan];
								w = gvx[pan] + gvx[aux] * q; gvx[pan] = w; gvx[aux] = w;
								w = gvy[pan] + gvy[aux] * q; gvy[pan] = w; gvy[aux] = w;
								w = gvz[pan] + gvz[aux] * q; gvz[pan] = w; gvz[aux] = w;
							}
						}
					}
				}
			}
		}

     underflow:	if (!glow) pcopy (adaptor, adapted);

		/* Calcoli principali. */

		if (!wait&&!glow) {
			tripper ();
			g_forces ();
		}

		/* Rotazione dei bordi delle bolle (per evidenziarle). */

		bbangle += 0.0025;
		if (bbangle>M_PI_4) bbangle -= M_PI_4;

	}
}
