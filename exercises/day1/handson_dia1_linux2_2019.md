---
output:
  pdf_document: default
  html_document: default
---
## Curso de Iniciación a la Secuenciación Masiva
BU-ISCIII

### Práctica 2 día 1: Software, Nextflow y Singularity

17-21 Junio 2019, 7a Edición, Programa Formación Continua, ISCIII


#### Descripción
En esta práctica vamos a aprender a comprobar la versión de software que utilizamos y ejecutar pipelines en Nextflow sobre containers de Singularity para asegurarnos de que nuestro análsis es reproducible.

#### Notas importantes
* USA EL TABULADOR PARA GUIARTE EN LA TERMINAL Y AUTOCOMPLETAR NOMBRE DE RUTAS, NOMBRE DE ARCHIVOS Y COMANDOS. (“El tabulador es tu mejor aliado”)

* USA LOS CURSORES DEL TECLADO PARA MOVERTE POR EL HISTORIAL DE COMANDOS EJECUTADOS (PODRÁS VOLVER A USARLOS SIN NECESIDAD DE VOLVERLOS A ESCRIBIR).

* NO ES ACONSEJABLE USAR ESPACIOS, TILDES NI CARACTERES ESPECIALES, COMO LA "Ñ", AL PONER NOMBRES A FICHERO O DIRECTORIOS.

#### Ejercicios
1. Comprueba qué instalación de git estás usando, dónde se están buscando instalaciones, como qué usuario estás trabajando, la ayuda del programa, la versión del programa, y el manual de usuario.
```bash
which git
echo $PATH
echo $USER
git -h
git --version
man git
```

2. Ejecuta Nextflow por primera vez
```bash
nextflow
```

3. Para ejecutar un pipeline con Nextflow, simplemente dale la señal 'run' seguida del PATH al script principal:
```bash
nextflow run /home/$USER/Documents/introduction_to_bioinformatics
```

4. Ejecutar Nextflow con una configuración personalizada:
```bash
nextflow -C /home/$USER/Documents/introduction_to_bioinformatics/nextflow.config \
run /home/$USER/Documents/introduction_to_bioinformatics/main.nf
```

5. No es necesario tener los scripts en tu ordenador, puedes ejecutarlos directamente desde GitHub sin tener que descargar o instalar nada:
```bash
nextflow run BU-ISCIII/introduction_to_bioinformatics
```

6. Ejecuta el Nextflow del curso para ver la ayuda:
```bash
nextflow run BU-ISCIII/introduction_to_bioinformatics --help
```

7. Finalmente ejecuta el pipeline sobre el container de Singularity, donde está instalado todo el software que necesitaremos:
```bash
nextflow run BU-ISCIII/introduction_to_bioinformatics -profile singularity
```
