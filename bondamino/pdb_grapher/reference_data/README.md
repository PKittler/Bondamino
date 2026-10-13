# Referenzwerte für Bindungslängen und Bindungswinkel

## Rückgrat: Engh & Huber (2001)

Engh, R. A. & Huber, R. (2001). Structure quality and target parameters.
*International Tables for Crystallography*, Vol. F, ch. 18.3, 382–392.

Die Werte stehen als `BACKBONE_BOND_REFERENCES` und `BACKBONE_ANGLE_REFERENCES` in
`utils.py`. Übernommen sind sie samt Zuordnungsregeln (Gly, Pro, cis-Pro, Residuum vor
bzw. nach der Peptidbindung) aus der Datei `EnghHuber_IntlTblsF_1999.txt` des Programms
Dangle (Richardson-Labor, Duke University):

- https://github.com/rlabduke/javadev/blob/7a10fd2f4be7fa2f5162765d35bb8e975e9bce5b/dangle/resource/dangle/EnghHuber_IntlTblsF_1999.txt
- Lizenz: BSD-3-Clause, siehe `LICENSE-javadev.txt`

## Seitenketten: CCP4 Monomer Library

`ccp4_sidechains.json` enthält die Bindungslängen (`value_dist`, `value_dist_esd`) und
Bindungswinkel (`value_angle`, `value_angle_esd`) der 20 Standard-Aminosäuren und
gängiger modifizierter Aminosäuren ohne Wasserstoffatome. Übernommen ist nur, woran
mindestens ein Atom außerhalb des Rückgrats (N, CA, C, O, OXT) beteiligt ist, weil die
Bibliothek die Aminosäuren als freie Moleküle beschreibt.

- https://github.com/MonomerLibrary/monomers, Commit `1c97ea23a7df06373ed01cd7c3156a1f8e9b27e3`
- Lizenz: LGPL-3.0, siehe `LICENSE-monomers.txt` und `GPL-3.0.txt`

Neu erzeugen: die `.cif`-Dateien der in `build_ccp4_sidechains.py` aufgeführten
Residuen aus dem Repository (Ordner `a/` bis `z/`) in einen Ordner laden, dann

    python -I build_ccp4_sidechains.py <ordner> ccp4_sidechains.json
