/**
 * iSPoC4iGWOAMu = internetable Standalone Proof of Concept for iGWOAMu
 * iGWOAMu = internetable Grid Wise Obscenlightment with Alealetry and Multunicast unit
 * 
 * ----- Concepts clés -----
 * iSekNet = infrastructure Sekurit Netswerk = {{Spievel},{Chalot}}
 * Spievel ("miroir" en romanche) = la passerelle permettant à Alice de passer du Clair à l'Obscur grâce à sa clé CAT
 *                                  et à en revenir
 * Chalot (abri de jardin équivalent du sac-à-main : personne n'y trouve jamais rien sauf le proporiétaire) = le stockage
 * CAT = Crumb Allocation Table : pas de carte, pas de trésor.
 * 
 * ----- Principe clé d'architecture -----
 * 
 * Les hashcodes sont une source raisonnable d'aléa non brownien invalidant l'heuristique.
 * Même si la base est régulière, comme la rotation d'une clé n-aire.
 * 
 * ----- GWOAM et Bloc Solution Logiciel iGWOAMu -----
 * 
 * iGWOAMu = internetable Grid Wise Obscenlightment with Aleametry and Mult-unicast unit
 * 
 * C <-nm-> M <-vp-> VP <-nm-> T <-vn-> VN <-nm-> F <-vq-> VQ=O <-nm-> S <-sntp-> E=TVSOS
 *          |         |        |         |        |         |          |             |
 * aléa-----(---------+--------(---------+--------(---------+----------(-------------+
 *          |                  |                  |                    |             
 * métrie---+------------------+------------------+--------------------+
 * 
 * C <-nm-> M <-vp-> VP <-nm-> T <-vn-> VN <-nm-> F <-vq-> VQ=O <-nm-> S <-cccp-> E=TVSOS
 *          |         |        |         |        |         |          |             |
 * aléa-----(---------+--------(---------+--------(---------+----------(-------------+
 *          |                  |                  |                    |             
 * métrie---+------------------+------------------+--------------------+
 * 
 * ({s},{@}) <-{sntp,cccp}-> {flux obscurs parce que données pures et non informations cryptées}
 * SNTP = SekNet Transfer Protocol = protocole avec AVEC connexion
 * CCCP = Crumbs Concurrent Channeling Protocol = protocole SANS connexion (datagrammes)
 * 
 * C = Clair --- M = {miette} --- T = {tranche} --- F = {fragment} --- O = Obscur
 * V{N, P, Q} = Vernamisé {Négatif, Positif, Quignon}
 * 
 * ian_ = indice absolu du nit _
 * i_ = indice de la {miette, tranche} ou du {fragment, flux}
 * irn_ = indice relatif à i_ du nit _
 * 
 * rk = rotation de K --- hrk = hahscode de la rotation de K --- shark[n,p] = somme des hrk [négatifs, positifs]
 * n_ = nombre de _
 * r_ = rang de _ (sauf rk)
 * t_ = nombre de nits dans _
 * 
 * ----- Calculs -----
 * |.| = taille ou valeur absolue de .
 * ---
 * tC = |C|
 * tK = |K|
 * ---
 * K
 * rk[i]            = K[i..|K|[ . K[0..i[
 * hrk[]            = rk[i].hashcode
 * shark            = +_i{|hrk[i]|}
 * sharkn           = +_i{|hrk[i]|, hrk[i]<0}
 * sharkp           = +_i{|hrk[i]|, 0<hrk[i]}
 * ---
 * nm               = |K|
 * np               = |{hrk[i], 0<hrk[i]}|
 * nt               = |K|
 * nn               = |{hrk[i], hrk[i]<0}|
 * nf               = |K|
 * nq               = |K|
 * ns               = |K|
 * ne               = |K|
 * ---
 * tm[0..nm[        = ((int)tC/nm)+((i<(tC%nm))?1:0)
 * ---
 * rp[0..np[        = |{j<i, 0<hrk[j]<=hrk[i]}|+|{j>i, 0<hrk[j]<hrk[i]}|
 * tp[1..np[        = {(int)((long)hrk[i]*((int)tC/np))/sharkp, hrk[i]>0}
 * tp[0]            = tC-+_[1..np[{tp[i]}
 * ---
 * tt[0..nt[        = ((int)tC/nt)+((i<(tC%nt))?1:0)
 * ---
 * rn[0..nn[        = |{j<i, 0<hrk[i]<=hrk[j]}|+|{j>i, 0<hrk[i]<hrk[j]}|
 * tn[1..nn[        = {(int)((long)hrk[i]*((int)tC/nn))/sahrkn, hrk[i]<0}
 * tn[0]            = tC-+_[1..nn[{tn[i]}
 * ---
 * rf[0..nf]        = |{j<i, hrk[i%nf]<=hrk[j%nf]}|+|{j>i, hrk[i%nf]<hrk[j%nf]}|
 * tf[0..nf]        = ((int)tC/nf)+((i<(tC%nf))?1:0)
 * ---
 * rq[0..nq[        = |{j<i, hrk[i]<=hrk[j]}|+|{j>i, hrk[i]<hrk[j]}|
 * tq[1..nq[        = |{(int)((long)hrk[i]*((int)tC/nn))/shark}|
 * tq[0]            = tC-+_[1..nq[{tn[i]}
 * ---
 * ts[0..ns[        = ((int)tC/ns)+((i<(tC%ns))?1:0)
 * ---
 * te[0..ne[        = ((int)tC/ne)+((i<(tC%ne))?1:0)
 * re[0..ne]        = |{j<i, hrk[i%nq]<=hrk[j%nq]}|+|{j>i, hrk[i%nq]<hrk[j%nq]}|
 * 
 * ----- Implémentation -----
 * C <-nm-> M <-vp-> VP <-nm-> T <-vn-> VN <-nm-> F <-vq-> VQ=O <-nm-> S <-{sntp,cccp}-> E=TVSOS
 * M ianm - émietté
 * P ianmv - p-morcellé
 * T ianmvm - tranché
 * N ianmvmv - n-morcellé
 * F ianmvmvm - fragmenté
 * Q ianmvmvmv - q-morcellé
 * S ianmvmvmvm - segmenté
 * E ianmvmvmvmv  - envoyé
 * 
 * Classes clés : AtomicReferenceArray, ArrayList, StringBilder
 */
package igwoamu;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

/**
 *
 * @author vSekNet
 */
public class IGWOAMu {
    
    final static boolean TRACER = false;
    
    /*
    * C = Clair --- M = {miette} --- T = {tranche} --- F = {fragment} --- O = Obscur
    * V{N, P, Q} = Vernamisé {Négatif, Positif, Quignon} */
    static enum eTransformation {M, VN, T, VP, F, VQ, O, E};
    static StringBuilder C = new StringBuilder("Bachi-bouzouk ! Moule à gaufres ! Australopithèque !");
    int tC;
    static StringBuilder K = new StringBuilder("Haddock");
    int tK;
    static StringBuilder O;
    static StringBuilder É;
    
    /* ian_ = indice absolu du nit _
    * i_ = indice de la {miette, tranche} ou du {fragment, flux}
    * irn_ = indice relatif à i_ du nit _ */
    static int ianc, ianm, ianp, iant, iann, ianf, ianq, ians, iano;
    static int[] im, it, if_, is;
    static int[] irnm, irnt, irnf, irns;
    
    /* rk = rotation de K --- hrk = hahscode de la rotation de K --- shark[n,p] = somme des hrk [négatifs, positifs]
    * n_ = nombre de _
    * r_ = rang de _ (sauf rk)
    * t_ = nombre de nits dans _*/
    String[]  rk;
    int[] hrk;
    long shark=0, sharkp=0, sharkn=0;
    int nm, np, nt, nn, nf, nq, ns, ne;
    int[] rm, rp, rt, rn, rf, rq, rs, re;
    int[] tm, tp, tt, tn, tf, tq, ts, te;

    IGWOAMu() {
        /* ----- Calculs -----
        * |.| = taille ou valeur absolue de .*/
        if (TRACER) System.out.printf("\n\n");
        
        /*
        * tC = |C|
        * tK = |K|*/
        tC = C.length();
        tK = K.length();
        if (TRACER) System.out.printf("C=\"%s\"\ntC=%d\nK=\"%s\"\ntK=%d\n\n", C, tC, K, tK);
        
        /* K
        * rk[i]            = K[i..|K|[ . K[0..i[
        * hrk[]            = rk[i].hashcode
        * shark            = +_i{|hrk[i]|}
        * sharkn           = +_i{|hrk[i]|, hrk[i]<0}
        * sharkp           = +_i{|hrk[i]|, 0<hrk[i]} */
        rk = new String[K.length()];
        hrk = new int[K.length()];
        for (int i=0; i<K.length(); i++) {
            rk[i] = K.substring(i).concat(K.substring(0,i));
            hrk[i] = rk[i].hashCode();
            shark += Math.abs(hrk[i]);
            if (hrk[i]<0) sharkn += hrk[i];
            if (hrk[i]>0) sharkp += hrk[i];
        }
        if (TRACER) System.out.printf("rk = %s\nhrk = %s\nshark_ = %d %d %d/%d\n\n",Arrays.toString(rk), Arrays.toString(hrk), sharkn, sharkp, shark, (sharkp-sharkn));

        /* nm=nf=nt=ns=ne     = |K|
        * nq               = |{hrk[i]}| = |K|
        * nn               = |{hrk[i], hrk[i]<0}|
        * np               = |{hrk[i], 0<hrk[i]}|*/
        nm = K.length();
        np = 0;
        nt = K.length();
        nn = 0;
        nf = K.length();
        nq = K.length();
        ns = K.length();
        ne = K.length();
        for (int i=0; i<nq; i++) {
            if (hrk[i]<0) nn++;
            if (hrk[i]>0) np++;
        }
        if (TRACER) System.out.printf("nm=%d, np=%d, nt=%d, nn=%d, nf=%d, nq=%d, ns=%d, ne=%d\n\n", nm, np, nt, nn, nf, nq, ns, ne);
        
        /* tm[0..nm[        = ((int)tC/nm)+((i<(tC%nm))?1:0)*/
        tm = new int[nm];
        for (int i=0; i<nm; i++)
            tm[i] = ((int)tC/nm)+((i<(tC%nm))?1:0);
        if (TRACER) System.out.printf("tm[]=%s\n\n",Arrays.toString(tm));
        
        /* rp[0..np[        = |{j<i, 0<hrk[j]<=hrk[i]}|+|{j>i, 0<hrk[j]<hrk[i]}|
        * tp[1..np[        = {(int)((long)hrk[i]*((int)tC/np))/sharkp, hrk[i]>0}
        * tp[0]            = tC-+_[1..np[{tp[i]}*/
        rp = new int[np];
        tp = new int[np];
        for (int i=0; i<np; i++) {
            rp[i] = 0;
            tp[i] = 0;
        }
        int ip=0;
        for (int i=0; i<nq; i++) {
            if (hrk[i]>0) {
                for (int j=0; j<nq; j++) {
                    if ((j<i)&&(hrk[j]<=hrk[i])&&(hrk[j]>0)) rp[ip]++;
                    if ((j>i)&&(hrk[j]< hrk[i])&&(hrk[j]>0)) rp[ip]++;
                }
                ip++;
            }
        }
        ip = 0;
        for (int i=0; i<nq; i++) {
            if (hrk[i]>0) {
                tp[ip] = (int) ((double)hrk[i]*tC/sharkp);
                ip++;
            }
        }
        tp[0] = tC;
        for (int i=1; i<np; i++)
            tp[0] -= tp[i];
        if (TRACER) System.out.printf("np=%d\nrp[]=%s\ntp[]=%s\n\n",np,Arrays.toString(rp),Arrays.toString(tp));
        
        /* tt[0..nt[        = ((int)tC/nt)+((i<(tC%nt))?1:0)*/
        tt = new int[nt];
        for (int i=0; i<nt; i++)
            tt[i] = ((int)tC/nt)+((i<(tC%nt))?1:0);
        if (TRACER) System.out.printf("tt[]=%s\n\n",Arrays.toString(tt));
        
        /* rn[0..nn[        = |{j<i, 0<hrk[i]<=hrk[j]}|+|{j>i, 0<hrk[i]<hrk[j]}|
        * tn[1..nn[        = {(int)((long)hrk[i]*((int)tC/nn))/sahrkn, hrk[i]<0}
        * tn[0]            = tC-+_[1..nn[{tn[i]}*/
        rn = new int[nn];
        tn = new int[nn];
        for (int i=0; i<nn; i++) {
            rn[i] = 0;
            tn[i] = 0;
        }
        int in=0;
        for (int i=0; i<nq; i++) {
            if (hrk[i]<0) {
                for (int j=0; j<nq; j++) {
                    if ((j<i)&&(hrk[j]<=hrk[i])) rn[in]++;
                    if ((j>i)&&(hrk[j]< hrk[i])) rn[in]++;
                }
                in++;
            }
        }
        in = 0;
        for (int i=0; i<nq; i++) {
            if (hrk[i]<0) {
                tn[in] = (int) ((double)hrk[i]*tC/sharkn);
                in++;
            }
        }
        tn[0] = tC;
        for (int i=1; i<nn; i++)
            tn[0] -= tn[i];
        if (TRACER) System.out.printf("nn=%d\nrn[]=%s\ntn[]=%s\n\n",nn,Arrays.toString(rn),Arrays.toString(tn));
        
        /* rf[0..nf]        = |{j<i, hrk[i%nf]<=hrk[j%nf]}|+|{j>i, hrk[i%nf]<hrk[j%nf]}|
        * tf[0..nf]        = ((int)tC/nf)+((i<(tC*nf))?1:0)*/
        rf = new int[nf];
        tf = new int[nf];
        for (int i=0; i<nf; i++) {
            rf[i] = 0;
            for (int j=0; j<i; j++)
                if (hrk[i%nf]<=hrk[j%nf]) rf[i]++;
            for (int j=i; j<nf; j++)
                if (hrk[i%nf]< hrk[j%nf]) rf[i]++;
        }
        for (int i=0; i<nf; i++)
            tf[i] = ((int)tC/nf)+((i<(tC%nf))?1:0);
        if (TRACER) System.out.printf("rf[]=%s\ntf[]=%s\n\n",Arrays.toString(rf),Arrays.toString(tf));
        
        
        /* rq[0..nq[        = |{j<i, hrk[i]<=hrk[j]}|+|{j>i, hrk[i]<hrk[j]}|
        * tq[1..nq[        = |{(int)((long)hrk[i]*((int)tC/nq))/shark}|
        * tq[0]            = tC-+_[1..nq[{tn[i]}*/
        rq = new int[nq];
        tq = new int[nq];
        for (int i=0; i<nq; i++) {
            rq[i] = 0;
            tq[i] = 0;
        }
        for (int i=0; i<nq; i++) {
            for (int j=0; j<nq; j++) {
                if ((j<i)&&(hrk[j]<=hrk[i])) rq[i]++;
                if ((j>i)&&(hrk[j]< hrk[i])) rq[i]++;
            }
        }
        for (int i=0; i<nq; i++) {
            tq[i] = (int) ((double)Math.abs(hrk[i])*tC/shark);
        }
        tq[0] = tC;
        for (int i=1; i<nq; i++)
            tq[0] -= tq[i];
        if (TRACER) System.out.printf("nq=%d\nrq[]=%s\ntq[]=%s\n\n",nq,Arrays.toString(rq),Arrays.toString(tq));
        
        /* ts[0..nq[        = ((int)tC/ns)+((i<(tC%ns))?1:0)*/
        ts = new int[ns];
        for (int i=0; i<ns; i++)
            ts[i] = ((int)tC/ns)+((i<(tC%ns))?1:0);
        if (TRACER) System.out.printf("ns=%d\nts[]=%s\n\n",ns,Arrays.toString(ts));
        
        /* te[0..nq[        = ((int)tC/ne)+((i<(tC%ne))?1:0)
        * re[0..ne]        = |{j<i, hrk[i%nq]<=hrk[j%nq]}|+|{j>i, hrk[i%nq]<hrk[j%nq]}| */
        re = new int[ne];
        te = new int[ne];
        for (int i=0; i<ne; i++) {
            re[i] = 0;
            te[i] = 0;
        }
        for (int i=0; i<ne; i++) {
            for (int j=0; j<ne; j++) {
                if ((j<i)&&(hrk[j%nq]<=hrk[i%nq])) re[i]++;
                if ((j>i)&&(hrk[j%nq]< hrk[i%nq])) re[i]++;
            }
        }
        for (int i=0; i<ne; i++)
            te[i] = ((int)tC/ne)+((i<(tC%ne))?1:0);
        if (TRACER) System.out.printf("ne=%d\nre[]=%s\nte[]=%s\n\n",ne,Arrays.toString(re),Arrays.toString(te));
        
        //=====================================================================
        for (int ianc=0; ianc<tC; ianc++) {
            if ((iano4ianc(ianc)<0)||(ianc4iano(ianc)>=tC)) System.exit(-ianc);
        }
        if (TRACER) System.out.println("0=<iano4ianc(ianc)<tC");
        for (int ianc=0; ianc<tC; ianc++) {
            if (ianc != ianc4iano(iano4ianc(ianc))) System.exit(-ianc);
        }
        if (TRACER) System.out.println("0=<iano4ianc(ianc)<tC");
        if (TRACER) System.out.println("ianc == ianc4iano(iano4ianc(ianc))");
        for (int ianc=0; ianc<tC; ianc++) {
            for (int janc=0; janc<tC; janc++) {
                if ((ianc != janc) && (ianc4iano(iano4ianc(ianc)) == ianc4iano(iano4ianc(janc)))) System.exit(((ianc==0)?-tC:-ianc));
            }
        }
        System.out.println("0=<iano4ianc(ianc)<tC");
        System.out.println("ianc == ianc4iano(iano4ianc(ianc))");
        System.out.println("(ianc !=janc) => (ianc4iano(ianc) !=ianc4iano(janc))");
        
        //=====================================================================
        O = new StringBuilder(C);
        for (int ianc=0; ianc<tC; ianc++)
            O.setCharAt(iano4ianc(ianc), C.charAt(ianc));
        É = new StringBuilder(C);
        for (int iano=0; iano<tC; iano++)
            É.setCharAt(ianc4iano(iano), O.charAt(iano));
        System.out.println("C=\""+C+"\"\nO=\""+O+"\nE=");
        for (int i=0; i<ne; i++) {
            int in0 = 0;
            for (int j=0; j<ne; j++) if (re[j]<re[i]) in0+=te[j];
            System.out.print("   ");
            for (int j=0; j<in0; j++) System.out.print(" ");
            System.out.println(O.substring(in0, in0+te[i]));
        }
        for (int i=0; i<ne; i++) {
            int in0 = 0;
            for (int j=0; j<ne; j++) if (re[j]<re[i]) in0+=te[j];
            System.out.println("E["+i+"]=\""+O.substring(in0, in0+te[i])+"\"");
        }
        System.out.println("É=\""+É+"\"");
    }
    
    private int iano4ianc(int ianc) {
        // M
        int ianm = ((ianc%nm)*((int)tC/nm))+(((ianc%nm)<(tC%nm))?(ianc%nm):(tC%nm)) + ((int)ianc/nm);
        // P
        int irnm = ianm;
        int ianmv;
        int inv=0;
        for (int i=0; i<np; i++) {
            for (int j=0; j<np; j++) {
                if (rp[j] <rp[i]) inv += tp[j];
            }
            for (int j=0; j<np; j++) {
                if (j<i) irnm -= tp[j];
            }
        }
        ianmv = inv + irnm;
        // T
        int ianmvm = ((ianmv%nm)*((int)tC/np))+(((ianmv%np)<(tC%np))?(ianmv%np):(tC%np)) + ((int)ianmv/np);
        // N
        int ianmvmv;
        int irnmvm = ianmvm;
        int invv=0;
        for (int i=0; i<nn; i++) {
            for (int j=0; j<nn; j++) {
                if (rn[j] <rn[i]) invv += tn[j];
            }
            for (int j=0; j<nn; j++) {
                if (j<i) irnmvm -= tn[j];
            }
        }
        ianmvmv = invv + irnmvm;
        // F
        int ianmvmvm = ((ianmvmv%nf)*((int)tC/nf))+(((ianmvmv%nf)<(tC%nf))?(ianmvmv%nf):(tC%nf)) + ((int)ianmvmv/nf);
        // Q
        int ianmvmvmv;
        int irnmvmvm = ianmvmvm;
        int invvv=0;
        for (int i=0; i<nq; i++) {
            for (int j=0; j<nq; j++) {
                if (rq[j] <rq[i]) invvv += tq[j];
            }
            for (int j=0; j<nq; j++) {
                if (j<i) irnmvmvm -= tq[j];
            }
        }
        ianmvmvmv = invvv + irnmvmvm;
        // S
        int ianmvmvmvm = ((ianc%nm)*((int)tC/nm))+(((ianc%nm)<(tC%nm))?(ianc%nm):(tC%nm)) + ((int)ianc/nm);
        // E
        
        if (TRACER) System.out.printf("(%d,%d,%d,%d,%d,%d,%d,%d,) %d->%d->%d->%d->%d->%d->%d->%d\n", nm, np, nt, nn, nf, nq, ns, ne, ianc,ianm,ianmv,ianmvm,ianmvmv,ianmvmvm,ianmvmvmv,ianmvmvmvm);
        return ianmvmvmvm;
    }
    private int ianc4iano(int iano) {
        int ianc = 0;
        while (iano4ianc(ianc) != iano) ianc++;
        return ianc; // C'est long, hein :)
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        InputStreamReader isr = new InputStreamReader(System.in);
        BufferedReader br = new BufferedReader(isr);
        System.out.printf("C=\"%s\"\nK=\"%s\"\n\nPour les changer entrer de nouvelles valeurs (ou deux retours chariot) :\n",C.toString(),K.toString());    
        String sC = br.readLine();
        String sK = br.readLine();
        if (!sC.isEmpty()) {
            C.delete(0, C.length());
            C.append(sC);
        }
        if (!sK.isEmpty()) {
            K.delete(0, K.length());
            K.append(sK);
        }
        System.out.printf("Travail avec\nC=\"%s\"\nK=\"%s\"\n", C.toString(), K.toString());
        
        IGWOAMu igwoamu = new IGWOAMu();
    }
}
