## Curso de Iniciación a la Secuenciación Masiva
BU-ISCIII

### Práctica 2 día 2: Manejo y gestión de ficheros

17-28 Mayo 2021, 8a Edición, Programa Formación Continua, ISCIII


#### Descripción
Al copiar o mover ficheros y/o directorios dentro del sistema o desde servidores de ficheros (FTP, samba, ...), tanto los permisos como propietario o grupo de esos ficheros, no siempre son los correctos para poder manipularlos.

El objetivo de la práctica es realizar desde la línea de comandos cambios en las características de los ficheros y/o directorios para que nuestro usuario del sistema pueda procesarlos y/o visualizarlos.

#### Notas importantes
* USA EL TABULADOR PARA GUIARTE EN LA TERMINAL Y AUTOCOMPLETAR NOMBRE DE RUTAS, NOMBRE DE ARCHIVOS Y COMANDOS. (“El tabulador es tu mejor aliado”)

* USA LOS CURSORES DEL TECLADO PARA MOVERTE POR EL HISTORIAL DE COMANDOS EJECUTADOS (PODRÁS VOLVER A USARLOS SIN NECESIDAD DE VOLVERLOS A ESCRIBIR).

* NO ES ACONSEJABLE USAR ESPACIOS, TILDES NI CARACTERES ESPECIALES, COMO LA "Ñ", AL PONER NOMBRES A FICHERO O DIRECTORIOS.

* LOS COMANDOS BÁSICOS QUE SIEMPRE DEBES RECORDAR: *pwd cd ls mkdir mv rm rmdir less nano*

#### Ejercicios
1) Acceder al directorio ‘Documentos’ de nuestro Home y dentro crear un directorio que se llame ‘practica_permisos’. Dentro de ‘practica_permisos’ crear otro directorio que se llame ‘copia_etc’. Copiar todos los ficheros que empiecen por ‘ho’ desde el directorio ‘/etc/’, a el directorio ‘copia_etc’ (usad el parámetro -v para ver el proceso de copia).

```bash
cd
cd Documentos
pwd
mkdir practica_permisos
cd practica_permisos
mkdir copia_etc
cd copia_etc
cp -v /etc/ho* ./ # (El punto “.” Indica que se copie en el directorio en el que te encuentras ubicado.)
```

2) Listar el contenido y visualizar los permisos del directorio ‘copia_etc’.

```bash
pwd
ll #o ls -l /home/alumno/Documentos/practica_permisos/copia_etc
```

3) Añadir permisos de escritura al grupo y resto de usuarios, en los ficheros que empiecen por ‘ho’. Listar el contenido y visualizar los permisos para ver los cambios.

```bash
pwd
chmod go+w ho* #o chmod 666 ho*
ll #o ls -l /home/alumno/Documentos/practica_permisos/copia_etc
```

4) Cambiar el grupo de los ficheros que empiece por ‘ho’ a ‘bioinfo’ (recuerda usar el comando ‘sudo’ para poder efectuar los cambios). Listar el contenido y visualizar para ver los cambios.

```bash
pwd
chgrp bioinfo ho*
ll #o ls -l /home/alumno/Documentos/practica_permisos/copia_etc
```

5) Cambiar el propietario de los ficheros que empiece por ‘ho’ a ‘root’ (recuerda usar el comando ‘sudo’ para poder efectuar los cambios). Listar el contenido y visualizar para ver los cambios.

```bash
pwd
sudo chown root ho*
ll #o ls -l /home/alumno/Documentos/practica_permisos/copia_etc
```

6) Quitar todos los permisos al ‘resto de usuarios’ (tercer bloque de permisos), de todos los ficheros dentro de ‘copia_etc’. Listar el contenido y visualizar los permisos para ver los cambios.

```bash
pwd
sudo chmod o-rwx * #o chmod 660 *
ll #o ls -l /home/Documentos/practica_permisos/copia_etc
```

7) Eliminar el directorio ‘practica_permisos’ y todo su contenido (usa el parámetro -f para forzar la eliminación).

```bash
pwd
cd ../../ #o cd ; cd Documentos
pwd
ls
rm -rf practica_permisos
ls
```

#### Nota final
* Podéis practicar estos ejercicios en cualquier ordenador, no solo en la máquian virtual del curso.

* Estos comandos son universales y funcionan en toda máquina que corra linux y similares, incluyendo macs y WSL.

* La manera más sencilla de practicarlos sin instalar nada es vía http://www.webminal.org/. En esta web puedes crearte un usuario de forma gratuita y abrir una terminal en una máquina remota, todo a través de vuestro explorador web. También contiene tutoriales complementarios que os pueden servir para afianzar lo aprendido hoy o repasar los comandos cuando tengáis necesidad de usarlos.

* *NOTA:* http://www.webminal.org/ ha sufrido recientemente un incendio en sus servidores y la plataforma está off-line temporalmente. Como alternativa podéis usar http://copy.sh/v86/?profile=archlinux, aunque esta carece de tutoriales.

```
Visita Webminal o copy.sh para practicar en casa: http://www.webminal.org/ o http://copy.sh/v86/?profile=archlinux
```
