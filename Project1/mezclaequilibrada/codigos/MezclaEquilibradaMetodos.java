package mezclaequilibrada.codigos;

import java.io.File;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.List;
import java.util.LinkedList;
import java.util.Queue;

/**
 * Clase que contiene los metodos que realizan mezcla equilibrada.
 * @author David Calderon Jimenez
 * @author Humberto Ignacio Hernandez Olvera
 * @author I\u00f1aky Ordiales Caballero
 * @version 1.9
 * @since 22/Octubre/2020
 */

public class MezclaEquilibradaMetodos {
  static int llamadaAp = 1;
  String dir1 = "mezclaequilibrada/porApellido_Iteracion#"+llamadaAp;
  String F0ap = dir1+"/F0.txt";
  String F1ap = dir1+"/F1.txt";
  String F2ap = dir1+"/F2.txt";
  String F1apiter = dir1+"/F1_Iteraciones.txt";
  String F2apiter = dir1+"/F2_Iteraciones.txt";
  String tmpap = dir1+"/tmp.txt";

  static int llamadaNom = 1;
  String dir2 = "mezclaequilibrada/porNombre_Iteracion#"+llamadaNom;
  String F0nom = dir2+"/F0.txt";
  String F1nom = dir2+"/F1.txt";
  String F2nom = dir2+"/F2.txt";
  String F1nomiter = dir2+"/F1_Iteraciones.txt";
  String F2nomiter = dir2+"/F2_Iteraciones.txt";
  String tmpnom = dir2+"/tmp.txt";

  String extension = ".txt";
  int iteracion = 0;

  /**
   * Este metodo sirve para crear los archivos que se utilizaran a lo largo de la ejecucion del metodo de mezcla equilibrada.
   * @param porApellido El criterio de ordenamiento, true si es por apellido, false si es por nombre.
   */
  public void creaArchivos(boolean porApellido) {

    File pAp = null;
    File pNom = null;
    // File f0 = null;
    File f1 = null;
    File f2 = null;
    File f1itera = null;
    File f2itera = null;
    File tmp = null;

    if (porApellido) {
      pAp = new File(dir1);
      pAp.mkdir();
      // f0 = new File(F0ap);
      f1 = new File(F1ap);
      f2 = new File(F2ap);
      f1itera = new File(F1apiter);
      f2itera = new File(F2apiter);
      tmp = new File(tmpap);
    } else {
      pNom = new File(dir2);
      pNom.mkdir();
      // f0 = new File(F0nom);
      f1 = new File(F1nom);
      f2 = new File(F2nom);
      f1itera = new File(F1nomiter);
      f2itera = new File(F2nomiter);
      tmp = new File(tmpnom);
    try {
      // f0.createNewFile();
      f1.createNewFile();
      f2.createNewFile();
      f1itera.createNewFile();
      f2itera.createNewFile();
      tmp.createNewFile();
    } catch(Exception e) {
      e.printStackTrace();
    }
  } // Fin método
}

  /**
   * Este metodo sirve para limpiar los archivos que se crearon en una ejecucion anterior.
   * @param porApellido El criterio de ordenamiento, true si es por apellido, false si es por nombre.
   */
  public void limpiaArchivos(boolean porApellido) {
    if (porApellido) {
      // File f0 = new File(F0ap);
      File f1 = new File(F1ap);
      File f2 = new File(F2ap);
      File f1itera = new File(F1apiter);
      File f2itera = new File(F2apiter);
      if (f1.exists() && f2.exists() && f1itera.exists() && f2itera.exists()) {
        try {
          // FileWriter f00 = new FileWriter(F0ap);
          FileWriter f11 = new FileWriter(F1ap);
          FileWriter f22 = new FileWriter(F2ap);
          File f11itera = new File(F1apiter);
          File f22itera = new File(F2apiter);
        } catch(Exception e) {
          e.printStackTrace();
        }
      } // Fin if si existen

    } else {
      // File f0 = new File(F0nom);
      File f1 = new File(F1nom);
      File f2 = new File(F2nom);
      File f1itera = new File(F1nomiter);
      File f2itera = new File(F2nomiter);
      if (f1.exists() && f2.exists() && f1itera.exists() && f2itera.exists()) {
        try {
          // FileWriter f00 = new FileWriter(F0nom);
          FileWriter f11 = new FileWriter(F1nom);
          FileWriter f22 = new FileWriter(F2nom);
          FileWriter f11itera = new FileWriter(F1nomiter);
          FileWriter f22itera = new FileWriter(F2nomiter);
        } catch(Exception e) {
          e.printStackTrace();
        }
      } // Fin if si existen
    } // Fin else
  } // Fin metodo

  /**
   * Este metodo realiza una copia del archivo original, con el proposito de reutilizar el original con otro metodo.
   * @deprecated No se usa mas en este proyecto.
   * @param path La ruta en donde se encuentra el archivo original.
   * @param porApellido El criterio de ordenamiento, true si es por apellido, false si es por nombre.
   * @return Una cadena que indica la ruta del archivo clonado.
   */
  @Deprecated public String clonOriginal(String path, boolean porApellido) {

    BufferedReader input = null;
    String path0 = "";

    if(porApellido) {
      path0 = F0ap;
      try {
        FileReader copyF = new FileReader(path);
        input = new BufferedReader(copyF);
        FileWriter f0 = new FileWriter(path0);
        BufferedWriter escribef0 = new BufferedWriter(f0);
        String linea;
        while((linea = input.readLine()) != null) {
          escribef0.write(linea);
          escribef0.newLine();
        }
        escribef0.close();
      } catch(Exception e) {
        e.printStackTrace();
      } finally {
          try {
            if (input == null) {
              input.close();
            }
          } catch(Exception e) {
            e.printStackTrace();
          }
      }
    } else {
      path0 = F0nom;
      try {
        FileReader copyF = new FileReader(path);
        input = new BufferedReader(copyF);
        FileWriter f0 = new FileWriter(path0);
        BufferedWriter escribef0 = new BufferedWriter(f0);
        String linea;
        while((linea = input.readLine()) != null) {
          escribef0.write(linea);
          escribef0.newLine();
        }
        escribef0.close();
      } catch(Exception e) {
        e.printStackTrace();
      } finally {
          try {
            if (input == null) {
              input.close();
            }
          } catch(Exception e) {
            e.printStackTrace();
          }
      } // Fin finally
    } // Fin else
    return path0;
  } // Fin metodo

  /**
   * Metodo para contar el numero de registros que contiene el archivo a ordenar.
   * @param path La ruta del archivo a ordenar (F0).
   * @return El numero de registros leidos.
   */
  public int comprueba(String path) {
    BufferedReader input = null;

    int contador = 0;

    try {
      FileReader file = new FileReader(path);
      input = new BufferedReader(file);
      String linea;
      while ((linea = input.readLine()) != null) {
        contador++;
      }
      return contador;
    } catch(Exception e) {
      e.printStackTrace();
    } finally {
        try {
          if (input == null) {
            input.close();
          }
        } catch(Exception e) {
          e.printStackTrace();
        }
    }
    return contador;
  }

  /**
   * Este metodo realiza el acomodo de los bloques de F0 a F1 y F2. Este metodo solo se llama al principio del programa.
   * @param path0 La ruta del archivo F0.
   * @param porApellido El criterio de ordenamiento, true si es por apellido, false si es por nombre.
   * @return Un arreglo de cadenas en donde se tiene las rutas de los archivos F0, F1 y F2 respectivamente.
   */
  public String[] ordenaI1(String path0, boolean porApellido){

    String path1;
    String path2;
    String path1itera;
    String path2itera;

    if (porApellido) {
      path1 = F1ap;
      path2 = F2ap;
      path1itera = F1apiter;
      path2itera = F2apiter;
    } else {
      path1 = F1nom;
      path2 = F2nom;
      path1itera = F1nomiter;
      path2itera = F2nomiter;
    }

    BufferedReader input = null;
    int alterna = 1;
    int comparacion = -1;
    iteracion++;

    try {
      FileWriter f1itera = new FileWriter(path1itera);
      BufferedWriter escf1itera = new BufferedWriter(f1itera);
      escf1itera.write("----------- Iteracion #"+iteracion);
      escf1itera.newLine();
      escf1itera.close();
      FileWriter f2itera = new FileWriter(path2itera);
      BufferedWriter escf2itera = new BufferedWriter(f2itera);
      escf2itera.write("----------- Iteracion #"+iteracion);
      escf2itera.newLine();
      escf2itera.close();
    } catch(Exception e) {
      e.printStackTrace();
    }

    try {
      FileReader file0 = new FileReader(path0);
      input = new BufferedReader(file0);
      String linea;
      List<Persona> personas = new LinkedList<>();
      while ((linea = input.readLine()) != null) {
        Persona per = new Persona();
        String[] addto = linea.split(",");
        per.setNombre(addto[0]);
        per.setApellido(addto[1]);
        per.setNumCuenta(Long.parseLong(addto[2].trim()));
        if (!personas.isEmpty()) {
          if (porApellido) {
            comparacion = personas.get(personas.size()-1).getApellido().compareTo(per.getApellido());
          } else {
            comparacion = personas.get(personas.size()-1).getNombre().compareTo(per.getNombre());
          }
        } // Fin if not empty
        if (personas.isEmpty()) {
          personas.add(per);
        } else if(comparacion < 0) { // comparamos el ultimo con el leido, si es menor, se debe seguir agregando a personas
          personas.add(per);
        } else {
          personas = escribeArchivoF1F2(personas, per, alterna, path1, path2, path1itera, path2itera);
          alterna++;
        } // Fin de impresion en archivos
    } // Fin while que lee
    escribeArchivoFinal(personas, alterna, path1, path2, path1itera, path2itera); // Con este imprimimos el último bloque leído
    personas.clear();
    } catch(Exception e) {
        System.out.println("Algo anda mal");
    } finally {
      try {
        if (input == null) {
          input.close();
        }
      } catch(Exception e) {
        System.out.println("Algo anda mal");
      }
    } // Termina el finally, que termina del try-catch

    String[] auxPaths = {path0, path1, path2, path1itera, path2itera};
    return auxPaths;
  } // Fin metodo Pasa de F0 a F1 y F2;

  /**
   * Este metodo es el que realiza el ordenamiento, juntando cada linea de F1 con F2, intercalandolos entre si y escribiendolo en F0. Este metodo se repite tantas veces como sea necesario.
   * @param path0 La ruta del archivo F0.
   * @param path1 La ruta del archivo F1.
   * @param path2 La ruta del archivo F2.
   * @param porApellido El criterio de ordenamiento, true si es por apellido, false si es por nombre.
   * @param path1itera La ruta del archivo F1_Iteraciones.
   * @param path2itera La ruta del archivo F2_Iteraciones.
   * @return Un arreglo de cadenas en donde se tiene las rutas de los archivos F0, F1 y F2 respectivamente.
   */
  public String[] ordenaI2(String path0, String path1, String path2,boolean porApellido, String path1itera, String path2itera) {

    BufferedReader inputf1 = null;
    BufferedReader inputf2 = null;

    try {
      FileWriter f0 = new FileWriter(path0);  // Para limpiar el archivo F0
      BufferedWriter escribef0 = new BufferedWriter(f0);
      escribef0.close();
      List<Persona> personas = new LinkedList<>();

      FileReader filef1 = new FileReader(path1);
      inputf1 = new BufferedReader(filef1);
      FileReader filef2 = new FileReader(path2);
      inputf2 = new BufferedReader(filef2);

      String linef1 = null;
      String linef2 = null;
      int escrituras = 0;

      while ((linef1 = inputf1.readLine()) != null && (linef2 = inputf2.readLine()) != null) {
        Queue<Persona> personasf1 = new LinkedList<Persona>();
        Queue<Persona> personasf2 = new LinkedList<Persona>();

        personasf1 = creaCola(linef1.split("/"));
        personasf2 = creaCola(linef2.split("/"));

        boolean flag = personasf1.size() <= personasf2.size();

        if (flag) {
          personas = procesaColas(personasf1, personasf2, porApellido);
        } else {
          personas = procesaColas(personasf2, personasf1, porApellido);
        }

        f0 = new FileWriter(path0, true);
        escribef0 = new BufferedWriter(f0);

        for (Persona per : personas) {
          escribef0.write(per.impArchivoAux());
          escrituras++;
        }
        escribef0.newLine();
        escribef0.close();
        personas.clear();

      } // Fin while

      if (linef1 != null) {
        Queue<Persona> personasf1 = new LinkedList<Persona>();
        Queue<Persona> personasf2 = new LinkedList<Persona>();
        personasf1 = creaCola(linef1.split("/"));
        personas = procesaColas(personasf2, personasf1, porApellido);

        f0 = new FileWriter(path0, true);
        escribef0 = new BufferedWriter(f0);
        for (Persona per : personas) {
          escribef0.write(per.impArchivoAux());
          escrituras++;
        }
        escribef0.newLine();
        personas.clear();
        escribef0.close();
      }
      inputf1.close();
      inputf2.close();
    } catch(Exception e) {
      System.out.println("Algo salio mal");
      e.printStackTrace();
    }

    String[] pathsAux  = {path0, path1, path2, path1itera, path2itera};
    return pathsAux;
  } // Fin metodo

  /**
   * Este metodo crea los bloques que se escriben en F1 y F2 desde F0.
   * @param path0 La ruta del archivo F0.
   * @param path1 La ruta del archivo F1.
   * @param path2 La ruta del archivo F2.
   * @param path1itera La ruta del archivo F1_Iteraciones.
   * @param path2itera La ruta del archivo F2_Iteraciones.
   * @return Un arreglo de cadenas en donde se tiene las rutas de los archivos F0, F1 y F2 respectivamente.
   */
  public String[] ordenaI3(String path0, String path1, String path2, String path1itera, String path2itera) {  // Lee de F0 para modificar F1 y F2, separa por lineas de F0 en bloques

    int parametro = 1;
    BufferedReader inputf0 = null;
    iteracion++;

    try {   // Se elimina el contenido de F1 y F2, puesto que está dentro de F0 ordenado por bloque.
      FileWriter f1 = new FileWriter(path1);
      BufferedWriter escribef1 = new BufferedWriter(f1);
      escribef1.close();
      FileWriter f2 = new FileWriter(path2);
      BufferedWriter escribef2 = new BufferedWriter(f2);
      escribef2.close();

      FileReader filef0 = new FileReader(path0);
      inputf0 = new BufferedReader(filef0);

      String linef0;

      FileWriter f1itera = new FileWriter(path1itera, true);
      BufferedWriter escf1itera = new BufferedWriter(f1itera);
      escf1itera.newLine();
      escf1itera.write("----------Iteracion #"+iteracion);
      escf1itera.newLine();
      escf1itera.close();
      FileWriter f2itera = new FileWriter(path2itera, true);
      BufferedWriter escf2itera = new BufferedWriter(f2itera);
      escf2itera.newLine();
      escf2itera.write("----------Iteracion #"+iteracion);
      escf2itera.newLine();
      escf2itera.close();

      while ((linef0 = inputf0.readLine()) != null) {
        if (parametro % 2 != 0) {
          f1 = new FileWriter(path1, true);
          escribef1 = new BufferedWriter(f1);
          escribef1.write(linef0);
          escribef1.write("\n");
          escribef1.close();
          f1itera = new FileWriter(path1itera, true);
          escf1itera = new BufferedWriter(f1itera);
          escf1itera.write(linef0);
          escf1itera.write("\n");
          escf1itera.close();
          parametro++;
        } else {
          f2 = new FileWriter(path2, true);
          escribef2 = new BufferedWriter(f2);
          escribef2.write(linef0);
          escribef2.write("\n");
          escribef2.close();
          f2itera = new FileWriter(path2itera, true);
          escf2itera = new BufferedWriter(f2itera);
          escf2itera.write(linef0);
          escf2itera.newLine();
          escf2itera.close();
          parametro++;
        } // Fin if-else
      } // Fin del while
      inputf0.close();
    } catch(Exception e) {
      System.out.println("Algo salio mal");
    }

    String[] pathsAux = {path0, path1, path2, path1itera, path2itera};
    return pathsAux;
  } // Fin del metodo

  /**
   * Este metodo sirve para regresar el archivo con los registros ordenados al formato en que estos fueron leidos del archivo original, utilizando un archivo tmp.
   * @param path0 La ruta del archivo F0.
   * @param porApellido El criterio de ordenamiento, true si es por apellido, false si es por nombre.
   */
  public void regresaAFormato(String path0, boolean porApellido) {

    BufferedReader input = null;
    String pathtmp;

    if (porApellido) {
      pathtmp = tmpap;
    } else {
      pathtmp = tmpnom;
    }

    File f = new File(pathtmp);

    try {
    FileReader filef0 = new FileReader(path0);
    input = new BufferedReader(filef0);
    FileWriter escribetmpf = new FileWriter(pathtmp);
    String linea;
    while((linea = input.readLine()) != null) { // Copia la unica linea de F0 a un tmp.txt
      // linea = sc.nextLine();
      BufferedWriter fAux = new BufferedWriter(escribetmpf);
      fAux.write(linea);
      fAux.close();
      }
      input.close();
    } catch(Exception e) {
      e.printStackTrace();
    } // Fin del try-catch

    try { // F0 ya fue copiado, se puede sobreescribir
      FileWriter f0 = new FileWriter(path0);  // Se ha limpiado F0
      f0.close();
      FileReader filetmp = new FileReader(pathtmp);
      input = new BufferedReader(filetmp);
      FileWriter escribef0 = new FileWriter(path0, true);
      String linea;
      while ((linea = input.readLine()) != null) {
        String[] personas = linea.split("/");
        for (String persona : personas) {
          String[] creaAdanEva = persona.split(",");
            Persona sujetoprueba = new Persona();
            sujetoprueba.setNombre(creaAdanEva[0]);
            sujetoprueba.setApellido(creaAdanEva[1]);
            sujetoprueba.setNumCuenta(Long.parseLong(creaAdanEva[2].trim()));
            escribef0 = new FileWriter(path0, true);
            escribef0.write(sujetoprueba.toString());
            escribef0.close();
          } // Fin foreach de cada sujeto
      } // fin while de lectura
      input.close();
      f.delete();
    } catch(Exception e) {
      e.printStackTrace();
    } // Fin del try-catch

  } // Fin del metodo

  /**
   * Crea la cola de personas, procesando la linea leida del archivo.
   * @param linea La linea del archivo que contiene los datos de las personas.
   * @return Una cola que contiene a todas las personas que se encontraban en esa linea del archivo.
   */
  public Queue<Persona> creaCola(String[] linea) {
    Queue<Persona> personasf = new LinkedList<Persona>();

    if (linea.length == 0) {
      return personasf;
    }
    // Creando lista de personas por linea del archivo
    for (String straux : linea) {
      String[] addto1 = straux.split(",");
      Persona per = new Persona();
      per.setNombre(addto1[0]);
      per.setApellido(addto1[1]);
      per.setNumCuenta(Long.parseLong(addto1[2].trim()));
      personasf.add(per);
    }
    return personasf;
  } // Fin metodo

  /**
   * Realiza el ordenamiento entre las colas con personas que fueron leidas de los archivos.
   * @param colaMenor La cola con la menor cantidad de personas.
   * @param colaMayor La cola con la mayor cantidad de personas.
   * @param porApellido El criterio de ordenamiento, true si es por apellido, false si es por nombre.
   * @return Una lista con las personas ordenadas segun el criterio.
   */
  public List<Persona> procesaColas(Queue<Persona> colaMenor, Queue<Persona> colaMayor, boolean porApellido) {
    List<Persona> personas = new LinkedList<>();  // Se crea una lista de personas para agregarlas de forma ordenada
    int sizeMayor = colaMayor.size();
    int sizeMenor = colaMenor.size();

    if (colaMenor.isEmpty()) {
      personas.addAll(colaMayor);
      return personas;
    }

    for (int i = 0; i < sizeMayor; i++) {
      if (i < sizeMenor) {
        personas = agregaPersona(personas, colaMenor.poll(), porApellido);
      } // Fin del if
      personas = agregaPersona(personas, colaMayor.poll(), porApellido);

    } // Fin del for
    return personas;
  }

  /**
   * Se agrega a una persona a la lista de personas, en la posicion correspondiente segun el criterio de ordenamiento.
   * Recursivamente se llama a si misma, con personas sin el primer elemento siempre, de tal forma que el caso base se cumpla, en donde se sacaron todos los menores al elemento que se quiere agregar y el primero de personas será mayor.
   * @param personas La lista de personas, puede o no estar vacia.
   * @param firstP La primera persona que se encuentra en la cola desde donde se desencolo.
   * @param porApellido El criterio de ordenamiento, true si es por apellido, false si es por nombre.
   * @return La lista de personas con la nueva persona ordenada.
   */
  public List<Persona> agregaPersona(List<Persona> personas, Persona firstP, boolean porApellido) {  // Agrega la primera persona de la cola de donde se llamó
    int comparacion = -1;
    if (personas.isEmpty()) {
      personas.add(firstP);
      return personas;
    }
    if (porApellido) {
      comparacion = personas.get(0).getApellido().compareTo(firstP.getApellido());
    } else {
      comparacion = personas.get(0).getNombre().compareTo(firstP.getNombre());
    }

    if (comparacion >= 0) { // Si el que se compara es menor que el primero de personas, se agrega al inicio
      personas.add(0, firstP);
      return personas;
    } else {
      List<Persona> primerP = new LinkedList<Persona>();  // Crea una lista de personas vacias para almacenar el primero de personas
      primerP.add(personas.remove(0));
      primerP.addAll(agregaPersona(personas, firstP, porApellido));
      /* Recursivamente se llama a si misma, con personas sin el primer elemento siempre, de tal forma que el caso base se cumpla, en donde se sacaron todos los menores al elemento que se quiere agregar y el primero de personas será mayor.
      */
      return primerP;
    }
  }

  /**
   * Para poder escribir en los archivos F1 o F2, alternando cada bloque que se tiene.
   * @param personas La lista de personas a escribir.
   * @param per La ultima persona que se tiene guardada en memoria y que se debe agregar a personas despues de escribirlas en el archivo.
   * @param alterna El criterio para decidir si se escriben en F1 o F2.
   * @param path1 La ruta del archivo F1.
   * @param path2 La ruta del archivo F2.
   * @param path1itera La ruta del archivo F1_Iteraciones.
   * @param path2itera La ruta del archivo F2_Iteraciones.
   * @return Una lista de personas, que contiene a la ultima persona leida pero que es menor al ultimo del bloque.
   */
  public List<Persona> escribeArchivoF1F2(List<Persona> personas, Persona per,int alterna,String path1,String path2, String path1itera, String path2itera) {
    if (alterna % 2 != 0) {
      try { // Para archivo F1
        FileWriter f1 = new FileWriter(path1, true);
        BufferedWriter escribe = new BufferedWriter(f1);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        f1.close();
        FileWriter f1itera = new FileWriter(path1itera, true);
        escribe = new BufferedWriter(f1itera);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        f1itera.close();
        personas.clear();
        personas.add(per);  // Se agrega la persona leída a personas, que es menor que el ultimo del bloque
      } catch (IOException e) {
        System.out.println("Algo anda mal");
      }
      return personas;
    } else {
      try { // Para archivo F2
        FileWriter f2 = new FileWriter(path2, true);
        BufferedWriter escribe = new BufferedWriter(f2);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        f2.close();
        FileWriter f2itera = new FileWriter(path2itera, true);
        escribe = new BufferedWriter(f2itera);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        f2itera.close();
        personas.clear();
        personas.add(per);  // Se agrega la persona leída a personas, que es menor que el ultimo del bloque
      } catch (IOException e) {
        System.out.println("Algo anda mal");
      } // Fin try-catch
      return personas;
    } // Fin else que imprime en archivo f2
  } // Fin del metodo

  /**
   * Este metodo sirve para escribir en el archivo, al ultimo bloque leido desde F0.
   * @param personas La lista de personas a escribir en el archivo.
   * @param alterna El criterio para decidir si escribe en F1 o F2.
   * @param path1 La ruta del archivo F1.
   * @param path2 La ruta del archivo F2.
   * @param path1itera La ruta del archivo F1_Iteraciones.
   * @param path2itera La ruta del archivo F2_Iteraciones.
   */
  public void escribeArchivoFinal(List<Persona> personas,int alterna,String path1,String path2, String path1itera, String path2itera) {
    if (alterna % 2 != 0) {
      try { // Para archivo F1
        FileWriter f1 = new FileWriter(path1, true);
        BufferedWriter escribe = new BufferedWriter(f1);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        f1.close();
        FileWriter f1itera = new FileWriter(path1itera, true);
        escribe = new BufferedWriter(f1itera);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        f1itera.close();
        personas.clear();
      } catch (IOException e) {
        System.out.println("Algo anda mal");
      }
    } else {
      try { // Para archivo F2
        FileWriter f2 = new FileWriter(path2, true);
        BufferedWriter escribe = new BufferedWriter(f2);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        f2.close();
        FileWriter f2itera = new FileWriter(path2itera, true);
        escribe = new BufferedWriter(f2itera);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        f2itera.close();
        personas.clear();
      } catch (IOException e) {
        System.out.println("Algo anda mal");
      } // Fin try-catch
    } // Fin else que imprime en archivo f2
  } // Fin del metodo

  /**
   * Este metodo sirve para borrar los archivos F1 y F2.
   * @param pathf1 La ruta del archivo F1.
   * @param pathf2 La ruta del archivo F2.
   */
  public void borraArchivos(String pathf1, String pathf2) {
    File f1 = new File(pathf1);
    File f2 = new File(pathf2);

    f1.delete();
    f2.delete();
  }

}
