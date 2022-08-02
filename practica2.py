import regex as re


if __name__ == '__main__':
    #Leemos el nombre del fichero
    fichero = input()
    #Expresion regular para que termine en DNA.txt
    patron = r'(.*)DNA\.txt'
    er_patron = re.compile(patron)
    result = er_patron.fullmatch(fichero)
    #Si no termina en DNA.txt, se informa y se vuelve a pedir el nombre
    #Si se pasa nombre vacio se termina.
    while not result:
        print('El nombre de fichero debe terminar en DNA.txt')
        fichero = input()
        if len(fichero) == 0:
            exit()
        result = er_patron.fullmatch(fichero)
    #Intentamos abrir el fichero
    #Si no se puede, se informa y se pide otro y se vuelve a intentar.
    #Si se pasa nombre vacio se termina.
    abierto = False
    while not abierto:
        try:
            archivo = open(fichero)
            abierto = True
        except FileNotFoundError:
            print('El fichero no existe o no se puede abrir')
            fichero = input()
            result = er_patron.fullmatch(fichero)
            while not result:
                if len(fichero) == 0:
                    exit()
                print('El nombre de fichero debe terminar en DNA.txt')
                fichero = input()
                result = er_patron.fullmatch(fichero)
    #ER para la cabecera del formato FASTA
    cabecera = r'((>(.*?)     )(\d+) nt( fragment)?\n)'
    #ER para las cadenas de A;G;C;T del formato FASTA
    DNA = r'^([A-Z]{10} ){0,4}[A-Z]{0,10}\n'
    er_cabecera = re.compile(cabecera)
    er_DNA = re.compile(DNA)
    #Contador para contar las letras que llevamos, para ver que coincide con la cabecera.
    contador = 0
    #Para almacenar el nombre de la cabecera
    nombre = ''
    #Para ir almacenando la cadena que se guarda en un mapa para la parte 2.
    parte2 = ''
    num_linea = 0
    #Creamos nombre del fichero de las proteinas usando grupos.
    new_nombre = er_patron.sub(r'\1Protein.txt',fichero)
    archivo_new = open(new_nombre, 'w+')
    #Imprimo el nombre del fichero de salida para indicar su creacion
    print(new_nombre)
    #Para almacenar la cadena de A,G,C,U (ya traducida la T) para luego recorrer en grupos de 3
    concatenacion = ''
    #Para solo traducir y trabajar con las cadenas que cumplen formato FASTA
    valido = 0
    #Para ir almacenando la traduccion a proteina.
    #Necesito esta nueva variable para almacenar la sustitucion (usando ER) pues asi es mucho mas sencillo
    #formatearla en el formato FASTA adecuadamte.
    nuevaCadena = ''
    #Para guardar el numero de la cabecera y comparar luego.
    numero = 0
    #Para crear la cabecera FASTA de la proteina.
    titulo = ''
    #Para contar las letras en la traduccion y ponerla en la cabecera de la proteina.
    contador2 = 0
    #Para saber si estamos en la primera linea de letras A,G,C,T
    primera = 0
    #Almacena fragment
    fragment = ''
    mapa = {'AGG':'R', 'AGA':'R', 'AGC':'S', 'AGU':'S', 'AAG':'K', 'AAA': 'K', 'AAC': 'N',
            'AAU': 'N', 'ACG': 'T', 'ACA': 'T', 'ACC': 'T', 'ACU': 'T', 'AUG': 'M', 'AUA': 'I',
            'AUC': 'I', 'AUU': 'I', 'CGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGU': 'R',
            'CAG': 'Q', 'CAA': 'Q', 'CAC': 'H', 'CAU': 'H', 'CCG': 'P', 'CCA': 'P',
            'CCC': 'P', 'CCU': 'P', 'CUG': 'L', 'CUA': 'L', 'CUC': 'L', 'CUU': 'L',
            'UGG': 'W', 'UGC': 'C', 'UGU': 'C', 'UAC': 'Y', 'UAU': 'Y', 'UCG': 'S',
            'UCA': 'S', 'UCC': 'S', 'UCU': 'S', 'UUG': 'L', 'UUA': 'L', 'UUC': 'F',
            'UUU': 'F', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GGU': 'G', 'GAG': 'E',
            'GAA': 'E', 'GAC': 'D', 'GAU': 'D', 'GCG': 'A', 'GCA': 'A', 'GCC': 'A',
            'GCU': 'A', 'GUG': 'V', 'GUA': 'V', 'GUC': 'V', 'GUU':'V', 'UAA': '', 'UGA': '', 'UAG': ''}
    #Mapa para almacenar las relaciones de la parte 2
    mapaParte2 = {}
    #ER para sustituir las T
    sust = r'T'
    #ER codon final al final de linea
    final = r'(.*UAA$|.*UGA$|.*UAG$)'
    #ER codon final en medio
    mal_final = r'UAA|UAG|UGA'
    er_malFinal = re.compile(mal_final)
    er_final = re.compile(final)
    er_sust = re.compile(sust)
    #ER grupos de 3
    codon = r'([A-Z]{3})'
    #ER inicio
    comienzo = r'^AUG'
    er_comienzo = re.compile(comienzo)
    er_codon = re.compile(codon)
    for linea in archivo:
        num_linea = num_linea + 1
        result = er_cabecera.fullmatch(linea)
        result2 = er_DNA.fullmatch(linea)
        #Si es una linea en blanco y hay algo que traducir, se traduce:
        #Como traducimos aquí, algunos errores que se muestren, se mostrarán con el nº de linea de la linea en blanco...
        if linea == '\n':
            valido = 0
            if concatenacion != '':
                result_final = er_final.fullmatch(concatenacion)
                #Si hay codon de final al final, comprobamos que el nº es correcto
                #Si no lo hay, error tambien
                if result_final:
                    if contador != numero:
                        concatenacion = ''
                        parte2 = ''
                        valido = 0
                        print('El nº de nucleotidos no coincide con el indicado: linea ', num_linea)
                        continue
                else:
                    concatenacion = ''
                    parte2 = ''
                    valido = 0
                    print('No se termina en alguno de los codones de fin: linea ', num_linea)
                    continue
                #Buscamos codones
                for r in er_codon.finditer(concatenacion):
                    entrada = concatenacion[r.start():r.end()].replace(' ', '')
                    entrada = entrada.replace('\n', '')
                    #Si no esta en las llaves, es que no cumple FASTA...
                    if entrada not in mapa.keys():
                        nuevaCadena = ''
                        concatenacion = ''
                        parte2 = ''
                        print('El fichero no cumple el formato FASTA: linea ', num_linea)
                        break
                    malFinal_result = er_malFinal.fullmatch(entrada)
                    #Si tiene codon de final en posicion no final, mal.
                    if malFinal_result and contador2 < numero/3 - 1:
                        print('Hay un codon de parada en una posicion distinta al final: linea ', num_linea)
                        nuevaCadena = ''
                        concatenacion = ''
                        parte2 = ''
                        break
                    #Traducimos el codon
                    nuevaCadena = nuevaCadena + er_codon.sub(mapa[entrada], concatenacion[r.start():r.end()], 1)
                    contador2 = contador2 + 1
                    #Formateamos para tener grupos de 10 y lineas de 50
                    if contador2%50 == 0 and mapa[entrada] != '':
                        nuevaCadena = nuevaCadena + '\n'
                    elif contador2%10 == 0 and mapa[entrada] != '':
                        nuevaCadena = nuevaCadena + ' '
                #Si no ha habido algun error, creamos la cadena de la parte 2 e imprimimos en el fichero de salida
                if concatenacion != '':
                    parte2 = parte2 + '\nProteina: ' + str(contador2-1) + ' aa\n\n' + nuevaCadena + '\n\n'
                    mapaParte2[nombre] = parte2
                    parte2 = ''
                    titulo = titulo + str(contador2-1) + ' aa' + fragment + '\n'
                    archivo_new.write(titulo)
                    nuevaCadena = nuevaCadena + '\n' + '\n' + '\n'
                    archivo_new.write(nuevaCadena)
                    nuevaCadena = ''
                    concatenacion = ''
        #Si es la cabecera, cogemos los datos que necesitamos: nº, fragment, nombre, etc
        #Si el nº no es multiplo de 3, mal.
        #Activamos el bit de primera linea y el de validez.
        elif result:
            valido = 1
            contador = 0
            contador2 = 0
            primera = 1
            titulo = result.group(2)
            nombre = result.group(3)
            parte2 = parte2 + '\nNombre: ' + nombre + '\n'
            if result.group(5):
                fragment = ' fragment'
            else:
                fragment = ''
            numero = int(result.group(4))
            parte2 = parte2 + '\nADN: ' + str(numero) + ' nt\n\n'
            if numero%3 != 0:
                print('El nº de aminoacidos debe ser multiplo de 3: linea ', num_linea)
                valido = 0
                primera = 0

        #Si es una linea de A,G,C y T y es valido (porque su cabecera anterior esta justo antes y era correcta), leemos
        #Hacemos la sustitucion de T a U.
        #Y si es la primera linea, comprobamos que se empieza por AUG. Si no, mal
        #Vamos creando la salida para la aprte 2 tambien.
        elif result2:
            if valido:
                parte2 = parte2 + linea
                nueva = er_sust.sub('U', linea)
                nueva = nueva.strip().replace(' ', '')
                nueva = nueva.replace('\n', '')
                concatenacion = concatenacion + nueva
                if primera:
                    result_comienzo = er_comienzo.match(concatenacion)
                    if not result_comienzo:
                        print('No empieza por AUG: linea ', num_linea)
                        valido = 0
                        concatenacion = ''
                        parte2 = ''
                    primera = 0
                contador = contador + len(nueva)
        #No cumple formato FASTA
        else:
            valido = 0
            print('El fichero no cumple el formato FASTA: linea ', num_linea)


    #Hay que cerrar los ficheros, especialmente en el de las proteinas
    #Si no lo cerramos y comparamos los ficheros antes de acabar la ejecucion, podria faltar informacion por volcar.
    archivo.close()
    archivo_new.close()
    # Leemos ER mientras no sea vacia.
    ER = input()
    while len(ER) > 0:
        er_ER = re.compile(ER)
        coincidencias = 0
        #Recorremos las llaves del mapa buscando coincidencias
        llaves = mapaParte2.keys()
        for llave in llaves:
            ER_result = er_ER.fullmatch(llave)
            #Si hay match, imprimimos la entrada del mapa y formateamos como se nos pide.
            if ER_result:
                coincidencias = coincidencias + 1
                salida = '=======================================================\n' + mapaParte2[llave]
                print(salida)
        print('=======================================================\n\n' + str(coincidencias) + ' coincidencias\n')
        ER = input()
