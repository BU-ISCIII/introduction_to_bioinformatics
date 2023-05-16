## Curso de Iniciación a la Secuenciación Masiva
BU-ISCIII

### Práctica 2: Software, Nextflow y Singularity

17-28 Mayo 2021, 8a Edición, Programa Formación Continua, ISCIII


#### Descripción
En esta práctica vamos a aprender a manejar environments de conda y a comprobar la versión de software que utilizamos para asegurarnos de que nuestro análsis es reproducible. También como instalar un programa en el env y finalmente vamos a descagar los materiales del curso usando git.

#### Notas importantes
* USA EL TABULADOR PARA GUIARTE EN LA TERMINAL Y AUTOCOMPLETAR NOMBRE DE RUTAS, NOMBRE DE ARCHIVOS Y COMANDOS. (“El tabulador es tu mejor aliado”)

* USA LOS CURSORES DEL TECLADO PARA MOVERTE POR EL HISTORIAL DE COMANDOS EJECUTADOS (PODRÁS VOLVER A USARLOS SIN NECESIDAD DE VOLVERLOS A ESCRIBIR).

* NO ES ACONSEJABLE USAR ESPACIOS, TILDES NI CARACTERES ESPECIALES, COMO LA "Ñ", AL PONER NOMBRES A FICHERO O DIRECTORIOS.

* LOS COMANDOS BÁSICOS QUE SIEMPRE DEBES RECORDAR: *pwd cd ls mkdir mv rm rmdir less nano*


#### Ejercicios


1. Ejecuta Conda por primera vez:
```bash
cd
conda
conda --version
which conda
```
¿Sabrías decir dónde está instalado el programa conda y su versión?

2. Investiguemos qué envs tenemos disponibles:
```bash
conda info --envs # o conda env list, es lo mismo
```
¿Estamos ahora mismo trabajando en algún env?

3. Cambiar de env:
```bash
conda activate ngs_course
conda env list
```
¿Cómo podemos saber que hemos cambiado de env?

4. Vamos a crear un env nuevo:
```bash
conda create --name git
```
¿En qué env estamos ahora?

5. Veamos qué programas tenemos instalados en el env:
```bash
conda list
```

6. Instalemos un programa en nuestro nuevo env:
```bash
conda install -c anaconda git
```
¿Cómo podemos ver que el nuevo programa ha sido instalado en el env?

7. Comprueba qué instalación de git estás usando, dónde se están buscando instalaciones, como qué usuario estás trabajando, la ayuda del programa, la versión del programa, y el manual de usuario.
```bash
which git
echo $PATH
echo $USER
git -h
git --version
man git
```
¿Cómo consultarías los manuales y ayudas de los comandos que ya conoces?

8. Usemos git para descargarnos los materials del curso:
```bash
git clone https://github.com/BU-ISCIII/introduction_to_bioinformatics.git
```
¿Cómo saldríais del env?
