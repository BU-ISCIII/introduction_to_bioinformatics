## Curso de Iniciación a la Secuenciación Masiva

BU-ISCIII

### Práctica 1: Comandos básicos de linux

#### Descripción

En esta práctica, se usaran los comandos más básicos para poder trabajar desde la terminal.

El objetivo es realizar desde la línea de comandos, aquellas tareas cotidianas que se realizan desde entorno grafico como crear, copiar, mover, visualizar o editar archivos.

#### Notas importantes

* USA EL TABULADOR PARA GUIARTE EN LA TERMINAL Y AUTOCOMPLETAR NOMBRE DE RUTAS, NOMBRE DE ARCHIVOS Y COMANDOS. (“El tabulador es tu mejor aliado”)

* USA LOS CURSORES DEL TECLADO PARA MOVERTE POR EL HISTORIAL DE COMANDOS EJECUTADOS (PODRÁS VOLVER A USARLOS SIN NECESIDAD DE VOLVERLOS A ESCRIBIR).

* NO ES ACONSEJABLE USAR ESPACIOS, TILDES NI CARACTERES ESPECIALES, COMO LA "Ñ", AL PONER NOMBRES A FICHERO O DIRECTORIOS.

* LOS COMANDOS BÁSICOS QUE SIEMPRE DEBES RECORDAR: *pwd cd ls mkdir mv rm rmdir less nano*

#### Cheats

| Diectories | Files |
|---|---|
| mkdir | touch |
| rmdir | ¿rm? |
| cd | less |
| | cat |
| | nano |

#### Ejercicios

1. Input y output en la terminal. Muestra el calendario actual, el calendario del mes de diciembre y la fecha y hora actual. Borra la pantalla con el comando clear o usando la combinación de teclas ctrl + l.

```bash
cal
cal -m 6 2021
date
clear
```

¿Se te ocurre como mostrar tu fecha de nacimiento? ¿Qué utilidad puede tener imprimir una fecha?

2. Prueba a usar los comandos who, whoami e id para ver la información que muestran.

```bash
who
whoami
id
uname -a
```

¿Sabes decir qué información muestra cada uno?

3. Abrir otro terminal y escribir sudo su para cambiar a usuario root (tendrás que poner la contraseña – OJO: La contraseña no se visualiza mientras la escribes) y vuelve a usar los comandos anteriores. Veras que la información ahora pertenece al usuario administrador. Salir del terminal en el que estés logueado como administrador (root). Para ello escribe exit (para volver al usuario alumno) y otra vez exit para cerrar la terminal.

```bash
sudo su
who
whoami
id
uname -a
exit
```

 ¿Qué ha cambiado respeco al hacer esto mismo con tu usario?

4. Abrir una nueva terminal y comprobar en que directorio estas ubicado con el comando pwd. Mostrar los archivos en esa carpeta usando el comando ls, y más información con los parámetros -l, -a, y -la.

```bash
pwd
ls
ls -l
ls -a
ls -la
```

¿Dónde estás? ¿Qué ficheros ves? ¿Y qué carpetas? ¿Qué hace el -a?

5. Desde tu home acceder (cd) al directorio ‘Documentos’ y dentro crear (mkdir) un directorio que se llame ‘practica_comandos’. Acceder dentro del directorio ‘practica_comandos’ y crear 2 directorios, uno se llamara ‘dir1’ y el otro ‘dir2’. Acceder dentro del directorio ‘dir1’ y crear un fichero de ‘texto’ vacio (touch o >) llamado ‘archivo1.txt’. Volver al directorio de inicio (/home/alumno).

```bash
cd Documentos
pwd
mkdir practica_comandos
cd practica_comandos
mkdir dir1 dir2
cd dir1 #o cd /home/alumno/Documentos/practica_comandos/dir1
> archivo1.txt #o touch archivo1.txt
cd #o cd /home/alumno/
```

¿Cuál sería el path relativo al archivo que acabas de crear? ¿Y el absoluto? ¿Cómo lo crearías desde tu home?

6. Usos del comando cat:

  Acceder (cd) al directorio ‘dir1’, crea un fichero ‘practica_cat.txt’ usando el editor nano (nano practica_cat.txt), añadir texto y guardar (pulsar ctrl + o), y sal del editor (ctrl + x). Listar (ls) el contenido de ‘dir1’ para comprobar que se ha creado. Visualizar en pantalla el contenido de ‘practica_cat.txt’ (cat practica_cat.txt).

  Hacer una copia (cp) de ‘practica_cat.txt’ dentro del directorio ‘dir1’ y cambia el nombre por ‘concatenar.txt’. Edita ‘concatenar.txt’ con nano y cambia algo su contenido. Listar el contenido de ‘dir1’ para comprobar que se ha creado. Visualizar en pantalla el contenido de concatenar.txt’. Ahora, visualizar en pantalla el contenido de ‘practica_cat.txt’ y ‘concatenar.txt’ con una sola instrucción.

  Por último, usar cat para guardar el contenido de los 2 ficheros en uno nuevo y llámalo ‘juntar_ficheros.txt’. Visualizar en pantalla. Limpiar pantalla (clear).

```bash
cd Documentos/practica_comandos/dir1
nano practica_cat.txt
# (ctrl + o) para guardar cambios
# (ctrl + x) para salir
ls
cat practica_cat.txt

cp practica_cat.txt concatenar.txt
nano concatenar.txt
# (ctrl + o) para guardar cambios
# (ctrl + x) para salir
ls
cat concatenar.txt

cat practica_cat.txt concatenar.txt #o cat *
cat practica_cat.txt concatenar.txt > juntar_ficheros.txt
cat -n juntar_ficheros.txt #(-n permite ver el número de líneas)
clear
```

Antes usamos > para algo diferente, ¿qué hace exactamente el comando >? ¿Dónde hemos creado el archivo concatenado? ¿Cómo podemos volver a home?

7. Comprobar que estas en el directorio ‘dir1’ con pwd y listar su contenido (ls). Listar el mismo directorio (´dir1´) pero esta vez para ver información detallada (permisos, propietario, grupo...) (ls -l).Prueba con otros parámetros del comando ls (-t, -S, -r).

```bash
pwd
ls
ls -l
ls -t
ls -S
ls -r
```

¿Lo que nos muestra pwd es ruta absoluta o relativa?

8. Desde ‘dir1’ listar el contenido del home del alumno (ls /home/alumno). Listar el mismo directorio pero esta vez para ver información de permisos (ls -l, o usar el alias ll que equivale a lo mismo) y por último vuelve a listarlo para ver permisos y archivos ocultos (ls -la). Recuerda que los archivos ocultos empiezan por un “.” y suelen ser archivos de configuración del sistema o preferencias de usuario. Prueba a usar rutas absolutas y rutas relativas.

```bash
pwd
ls /home/alumno #o ls ../../../../alumno/
ls -l /home/alumno #o ll ../../../../alumno/
ls -la /home/alumno #o ll -a ../../../../alumno
```

¿Para qué es interesante crear o ver los archivos ocultos?

9. Prueba a moverte (cd) y listar (ls) otros directorios del árbol de ficheros de Linux. Recuerda que para saber en qué directorio te encuentras se usa el comando pwd y que con el comando cd (sin parámetros) vuelves a tu directorio personal (/home/alumno).

```bash
cd /
pwd
ls

cd etc
pwd
ls

cd
pwd
ls
```

¿Cómo puedes ver todas las opciones para cambiar de destino?

10. Acceder al directorio ‘dir2’ (anteriormente creado) y mueve el contenido de ‘dir1’ a ‘dir2’ (mv). Listar ambos directorios por separado para comprobar que se han movido los ficheros. Cambiar el nombre al fichero ‘archivo1.txt’ por ‘archivo_importante.txt’ (mv). Visualiza con cat el fichero ‘passwd’ que se encuentra en la carpeta del sistema /etc y redireccionar la salida del fichero a otro fichero que se llame ‘archivo_importante.txt’ (cat /etc/passwd > archivo_importante.txt) (El fichero ‘passwd’ contiene información de los usuarios del sistema). Visualizar el contenido usando cat y less (recuerda que less paginan el contenido y que para salir se usa la tecla q [quit]).

```bash
cd /home/alumno/Documentos/practica_comandos/dir2
mv ../dir1/* ./ #o mv /home/alumno/Documentos/practica_comandos/dir1/* /home/alumno/Documentos/practica_comandos/dir2/
# Recuerda que el punto “.” se usa para indicar el directorio actual y el asterisco * es un carácter especial que indica cualquier cadena de caracteres y en este caso indica todo el contenido de dir1.
# OJO! si no pones el * o indicas el fichero, SE MOVERA EL DIRECTORIO ENTERO.
ls .
ls ../dir1/

mv archivo1.txt archivo_importante.txt
cat /etc/passwd > archivo_importante.txt
less archivo_importante.txt
# q para salir
```

¿Por qué usar less cuando tenemos word?

11. Ahora copiar (cp) el fichero ‘passwd’ dentro de ‘dir2’ con el nombre ‘listado_usuarios.txt’. Listar el directorio para ver los permisos y propietario de los ficheros (ls -l). Mostrar el contenido de ambos ficheros. Para salir de la visualización de less pulsa "q".

```bash
cp /etc/passwd listado_usuarios.txt #o cp /etc/passwd ./; mv passwd listado_usuarios.txt
ls -l #o ls -l /home/alumno/Documentos/practica_comandos/dir2
less listado_usuarios.txt # Pulsa q para salir
# Acuérdate que puedes limpiar la pantalla con el comando clear o pulsando ctrl + l
```

¿Para qué sirve less si ya sé usar cat?

12. Desde dentro del directorio ‘practica_comandos’ hacer una copia del directorio ‘dir2’ con nombre ‘copia_dir2’ (cp –r).

  Desde ‘practica_comandos’, eliminar (rm -r) todos los archivos (en este caso los directorios ‘dir1’ y ‘dir2’) cuyo nombre empiece por la letra d. Listar para comprobar que se ha eliminado.

  Eliminar todo el contenido del directorio ‘practica_comandos’ (rm –r).

  Por último, eliminar el directorio ‘practica_comandos’ usando comando rmdir (solo funciona si el directorio esta vacío). NOTA: usar el comando history para ver todos los comandos usados en la práctica.

```bash
cd ..
pwd #o si no estás en el directorio practica_comandos ir a él (cd)
ls
cp -r dir2 copia_dir2
ls
ls dir2
ls copia_dir2

rm -rv d*
ls

rm -r *
ls

cd ..
pwd
ls

rmdir practica_comandos
ls
```

#### Nota final

* Podéis practicar estos ejercicios en cualquier ordenador, no solo en la máquian virtual del curso.

* Estos comandos son universales y funcionan en toda máquina que corra linux y similares, incluyendo macs y WSL.

* La manera más sencilla de practicarlos sin instalar nada es vía <http://www.webminal.org/>. En esta web puedes crearte un usuario de forma gratuita y abrir una terminal en una máquina remota, todo a través de vuestro explorador web. También contiene tutoriales complementarios que os pueden servir para afianzar lo aprendido hoy o repasar los comandos cuando tengáis necesidad de usarlos.

* *NOTA:* <http://www.webminal.org/> ha sufrido recientemente un incendio en sus servidores y la plataforma está off-line temporalmente. Como alternativa podéis usar <http://copy.sh/v86/?profile=archlinux>, aunque esta carece de tutoriales.

```
Visita Webminal o copy.sh para practicar en casa: http://www.webminal.org/ o http://copy.sh/v86/?profile=archlinux
```
