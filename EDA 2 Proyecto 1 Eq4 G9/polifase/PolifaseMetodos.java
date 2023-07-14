package polifase;

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
 * Clase que contiene los metodos que realizan polifase.
 * @author David Calderon Jimenez
 * @author Humberto Ignacio Hernandez Olvera
 * @author I\u00f1aky Ordiales Caballero
 * @version 1.5
 * @since 29/Octubre/2020
 */

public class PolifaseMetodos {
  String dir1 = "porApellido";
  String dir2 = "porNombre";
  String g1 = "g1";
  String g2 = "g2";
  String tmp = "tmp";
  String extension = ".txt";
  int iteracion = 0;

  /**
   * Este metodo sirve para crear los archivos que se utilizaran a lo largo de la ejecucion del metodo de polifase.
   * @param porApellido El criterio de ordenamiento, true si es por apellido, false si es por nombre.
   */
  public void creaArchivos(boolean porApellido) {

    File pAp = null;
    File pNom = null;
    File g1 = null;
    File g2 = null;
    File g1itera = null;
    File g2itera = null;
    File tmp = null;

    if (porApellido) {
      pAp = new File("polifase/iteraciones/porApellido");
      pAp.mkdir();
      g1 = new File("polifase/iteraciones/porApellido/g1.txt");
      g2 = new File("polifase/iteraciones/porApellido/g2.txt");
      g1itera = new File("polifase/iteraciones/porApellido/g1_Iteraciones.txt");
      g2itera = new File("polifase/iteraciones/porApellido/g2_Iteraciones.txt");
      tmp = new File("polifase/iteraciones/porApellido/tmp.txt");
    } else {
      pNom = new File("polifase/iteraciones/porNombre");
      pNom.mkdir();
      g1 = new File("polifase/iteraciones/porNombre/g1.txt");
      g2 = new File("polifase/iteraciones/porNombre/g2.txt");
      g1itera = new File("polifase/iteraciones/porNombre/g1_Iteraciones.txt");
      g2itera = new File("polifase/iteraciones/porNombre/g2_Iteraciones.txt");
      tmp = new File("polifase/iteraciones/porNombre/tmp.txt");
    try {
      g1.createNewFile();
      g2.createNewFile();
      g1itera.createNewFile();
      g2itera.createNewFile();
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
  public void borraArchivos(boolean porApellido) {
    if (porApellido) {
      File g1 = new File("polifase/iteraciones/porApellido/g1.txt");
      File g2 = new File("polifase/iteraciones/porApellido/g2.txt");
      File g1itera = new File("polifase/iteraciones/porApellido/g1_Iteraciones.txt");
      File g2itera = new File("polifase/iteraciones/porApellido/g2_Iteraciones.txt");
      if (g1.exists() && g2.exists() && g1itera.exists() && g2itera.exists()) {
        try {
          FileWriter g11 = new FileWriter("polifase/iteraciones/porApellido/g1.txt");
          FileWriter g22 = new FileWriter("polifase/iteraciones/porApellido/g2.txt");
          File g11itera = new File("polifase/iteraciones/porApellido/g1_Iteraciones.txt");
          File g22itera = new File("polifase/iteraciones/porApellido/g2_Iteraciones.txt");
        } catch(Exception e) {
          e.printStackTrace();
        }
      } // Fin if si existen

    } else {
      File g1 = new File("polifase/iteraciones/porNombre/g1.txt");
      File g2 = new File("polifase/iteraciones/porNombre/g2.txt");
      File g1itera = new File("polifase/iteraciones/porNombre/g1_Iteraciones.txt");
      File g2itera = new File("polifase/iteraciones/porNombre/g2_Iteraciones.txt");
      if (g1.exists() && g2.exists() && g1itera.exists() && g2itera.exists()) {
        try {
          FileWriter g11 = new FileWriter("polifase/iteraciones/porNombre/g1.txt");
          FileWriter g22 = new FileWriter("polifase/iteraciones/porNombre/g2.txt");
          FileWriter g11itera = new FileWriter("polifase/iteraciones/porNombre/g1_Iteraciones.txt");
          FileWriter g22itera = new FileWriter("polifase/iteraciones/porNombre/g2_Iteraciones.txt");
        } catch(Exception e) {
          e.printStackTrace();
        }
      } // Fin if si existen
    } // Fin else
  } // Fin metodo

  /**
   * Metodo para contar el numero de registros que contiene el archivo a ordenar.
   * @param path La ruta del archivo a ordenar (dd).
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
   * Este metodo realiza el el ordenamiento interno con bloques de tamaño 4. Este metodo solo se llama al principio del programa.
   * @param path0 La ruta del archivo dd.
   * @param porApellido El criterio de ordenamiento, true si es por apellido, false si es por nombre.
   * @return Un arreglo de cadenas en donde se tiene las rutas de los archivos dd, g1 y g2 respectivamente.
   */
  public String[] ordenaI1(String path0, boolean porApellido){

    String path1;
    String path2;
    String path1itera;
    String path2itera;

    if (porApellido) {
      path1 = "polifase/iteraciones/porApellido/g1.txt";
      path2 = "polifase/iteraciones/porApellido/g2.txt";
      path1itera = "polifase/iteraciones/porApellido/g1_Iteraciones.txt";
      path2itera = "polifase/iteraciones/porApellido/g2_Iteraciones.txt";
    } else {
      path1 = "polifase/iteraciones/porNombre/g1.txt";
      path2 = "polifase/iteraciones/porNombre/g2.txt";
      path1itera = "polifase/iteraciones/porNombre/g1_Iteraciones.txt";
      path2itera = "polifase/iteraciones/porNombre/g2_Iteraciones.txt";
    }

    BufferedReader input = null;
    int alterna = 1;
    int comparacion = 1;
    iteracion++;

    try {
      FileWriter g1itera = new FileWriter(path1itera);
      BufferedWriter escg1itera = new BufferedWriter(g1itera);
      escg1itera.write("----------- Iteracion #"+iteracion);
      escg1itera.newLine();
      escg1itera.close();
      FileWriter g2itera = new FileWriter(path2itera);
      BufferedWriter escg2itera = new BufferedWriter(g2itera);
      escg2itera.write("----------- Iteracion #"+iteracion);
      escg2itera.newLine();
      escg2itera.close();
    } catch(Exception e) {
      e.printStackTrace();
    }

    try {
      FileReader file0 = new FileReader(path0);
      input = new BufferedReader(file0);
      String linea;
      List<Persona> personas = new LinkedList<>();
      int politera =  0;
      while ((linea = input.readLine()) != null) {
        
        int pos = 0;
        Persona per = new Persona();
        Persona sw = new Persona();
        String[] addto = linea.split(",");
        per.setNombre(addto[0]);
        per.setApellido(addto[1]);
        per.setNumCuenta(Long.parseLong(addto[2].trim()));

        if (porApellido) {
          //ALGORITMO DE ORDENAMIENTO INTERNO BUBBLE SORT
            int n = personas.size(); 
            for (int i = 0; i < n-1; i++){
              for (int j = 0; j < n-i-1; j++){
                if (0 <= personas.get(j).getApellido().compareTo(personas.get(j+1).getApellido())){
             
                    sw = personas.get(j); 
                    personas.set(j,personas.get(j+1)); 
                    personas.set(j+1,sw); 
                } 
              }
            }
        
          } else {
            int n = personas.size(); 
            for (int i = 0; i < n-1; i++){ 
              for (int j = 0; j < n-i-1; j++){
                if (0 <= personas.get(j).getNombre().compareTo(personas.get(j+1).getNombre())){
             
                    sw = personas.get(j); 
                    personas.set(j,personas.get(j+1)); 
                    personas.set(j+1,sw); 
                }
              }
            }
          }
            
          personas.add(per);
          politera++;
        //} 
        //En 4 se refiere  tamaño de los bloque que se usaran para la polifase
        if(politera == 4) {
          
          personas = escribeArchivog1g2(personas, alterna, path1, path2, path1itera, path2itera);
          alterna++;
          politera = 0;
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
  } // Fin metodo Pasa de dd a g1 y g2;

  /**
   * Este metodo es el que realiza el ordenamiento, juntando cada linea de g1 con g2, intercalandolos entre si y escribiendolo en dd. Este metodo se repite tantas veces como sea necesario.
   * @param path0 La ruta del archivo dd.
   * @param path1 La ruta del archivo g1.
   * @param path2 La ruta del archivo g2.
   * @param porApellido El criterio de ordenamiento, true si es por apellido, false si es por nombre.
   * @param path1itera La ruta del archivo g1_Iteraciones.
   * @param path2itera La ruta del archivo g2_Iteraciones.
   * @return Un arreglo de cadenas en donde se tiene las rutas de los archivos dd, g1 y g2 respectivamente.
   */
  public String[] ordenaI2(String path0, String path1, String path2,boolean porApellido, String path1itera, String path2itera) {

    BufferedReader inputg1 = null;
    BufferedReader inputg2 = null;

    try {
      FileWriter dd = new FileWriter(path0);  // Para limpiar el archivo dd
      BufferedWriter escribedd = new BufferedWriter(dd);
      escribedd.close();
      List<Persona> personas = new LinkedList<>();

      FileReader fileg1 = new FileReader(path1);
      inputg1 = new BufferedReader(fileg1);
      FileReader fileg2 = new FileReader(path2);
      inputg2 = new BufferedReader(fileg2);

      String lineg1 = null;
      String lineg2 = null;
      int escrituras = 0;

      while ((lineg1 = inputg1.readLine()) != null && (lineg2 = inputg2.readLine()) != null) {
        Queue<Persona> personasg1 = new LinkedList<Persona>();
        Queue<Persona> personasg2 = new LinkedList<Persona>();

        personasg1 = creaCola(lineg1.split("/"));
        personasg2 = creaCola(lineg2.split("/"));

        boolean flag = personasg1.size() <= personasg2.size();

        if (flag) {
          personas = procesaColas(personasg1, personasg2, porApellido);
        } else {
          personas = procesaColas(personasg2, personasg1, porApellido);
        }

        dd = new FileWriter(path0, true);
        escribedd = new BufferedWriter(dd);

        for (Persona per : personas) {
          escribedd.write(per.impArchivoAux());
          escrituras++;
        }
        escribedd.newLine();
        escribedd.close();
        personas.clear();

      } // Fin while

      if (lineg1 != null) {
        Queue<Persona> personasg1 = new LinkedList<Persona>();
        Queue<Persona> personasg2 = new LinkedList<Persona>();
        personasg1 = creaCola(lineg1.split("/"));
        personas = procesaColas(personasg2, personasg1, porApellido);

        dd = new FileWriter(path0, true);
        escribedd = new BufferedWriter(dd);
        for (Persona per : personas) {
          escribedd.write(per.impArchivoAux());
          escrituras++;
        }
        escribedd.newLine();
        personas.clear();
        escribedd.close();
      }
        inputg1.close();
        inputg2.close();
    } catch(Exception e) {
      System.out.println("Exepciones porque el numero de datos es dividido exactamente entre 4 por lo cual no quedan datos extras");
      e.printStackTrace();

    }
    String[] pathsAux  = {path0, path1, path2, path1itera, path2itera};
    return pathsAux;
  } // Fin metodo

  /**
   * Este metodo crea los bloques que se escriben en g1 y g2 desde dd.
   * @param path0 La ruta del archivo dd.
   * @param path1 La ruta del archivo g1.
   * @param path2 La ruta del archivo g2.
   * @param path1itera La ruta del archivo g1_Iteraciones.
   * @param path2itera La ruta del archivo g2_Iteraciones.
   * @return Un arreglo de cadenas en donde se tiene las rutas de los archivos dd, g1 y g2 respectivamente.
   */
  public String[] ordenaI3(String path0, String path1, String path2, String path1itera, String path2itera) {  // Lee de dd para modificar g1 y g2, separa por lineas de dd en bloques

    int parametro = 1;
    BufferedReader inputdd = null;
    iteracion++;

    try {   // Se elimina el contenido de g1 y g2, puesto que está dentro de dd ordenado por bloque.
      FileWriter g1 = new FileWriter(path1);
      BufferedWriter escribeg1 = new BufferedWriter(g1);
      escribeg1.close();
      FileWriter g2 = new FileWriter(path2);
      BufferedWriter escribeg2 = new BufferedWriter(g2);
      escribeg2.close();

      FileReader filedd = new FileReader(path0);
      inputdd = new BufferedReader(filedd);

      String linedd;

      FileWriter g1itera = new FileWriter(path1itera, true);
      BufferedWriter escg1itera = new BufferedWriter(g1itera);
      escg1itera.newLine();
      escg1itera.write("----------Iteracion #"+iteracion);
      escg1itera.newLine();
      escg1itera.close();
      FileWriter g2itera = new FileWriter(path2itera, true);
      BufferedWriter escg2itera = new BufferedWriter(g2itera);
      escg2itera.newLine();
      escg2itera.write("----------Iteracion #"+iteracion);
      escg2itera.newLine();
      escg2itera.close();

      while ((linedd = inputdd.readLine()) != null) {
        if (parametro % 2 != 0) {
          g1 = new FileWriter(path1, true);
          escribeg1 = new BufferedWriter(g1);
          escribeg1.write(linedd);
          escribeg1.write("\n");
          escribeg1.close();
          g1itera = new FileWriter(path1itera, true);
          escg1itera = new BufferedWriter(g1itera);
          escg1itera.write(linedd);
          escg1itera.write("\n");
          escg1itera.close();
          parametro++;
        } else {
          g2 = new FileWriter(path2, true);
          escribeg2 = new BufferedWriter(g2);
          escribeg2.write(linedd);
          escribeg2.write("\n");
          escribeg2.close();
          g2itera = new FileWriter(path2itera, true);
          escg2itera = new BufferedWriter(g2itera);
          escg2itera.write(linedd);
          escg2itera.newLine();
          escg2itera.close();
          parametro++;
        } // Fin if-else
      } // Fin del while
    } catch(Exception e) {
      System.out.println("Algo salio mal");
    } finally {
      try {
        if (inputdd == null) {
          inputdd.close();
        }
      } catch(Exception e) {
        System.out.println("Algo anda mal");
      }
    } // Termina el finally, que termina del try-catch

    String[] pathsAux = {path0, path1, path2, path1itera, path2itera};
    return pathsAux;
  } // Fin del metodo

  /**
   * Este metodo sirve para regresar el archivo con los registros ordenados al formato en que estos fueron leidos del archivo original, utilizando un archivo tmp.
   * @param path0 La ruta del archivo dd.
   * @param porApellido El criterio de ordenamiento, true si es por apellido, false si es por nombre.
   */
  public void regresaAFormato(String path0, boolean porApellido) {

    BufferedReader input = null;
    String pathtmp;

    if (porApellido) {
      pathtmp = "polifase/iteraciones/porApellido/tmp.txt";
    } else {
      pathtmp = "polifase/iteraciones/porNombre/tmp.txt";
    }

    File f = new File(pathtmp);

    try {
    FileReader filedd = new FileReader(path0);
    input = new BufferedReader(filedd);
    FileWriter escribetmpf = new FileWriter(pathtmp);
    String linea;
    while((linea = input.readLine()) != null) { // Copia la unica linea de dd a un tmp.txt
      // linea = sc.nextLine();
      BufferedWriter fAux = new BufferedWriter(escribetmpf);
      fAux.write(linea);
      fAux.close();
      }
      input.close();
    } catch(Exception e) {
      e.printStackTrace();
    } // Fin del try-catch

    try { // dd ya fue copiado, se puede sobreescribir
      FileWriter dd = new FileWriter(path0);  // Se ha limpiado dd
      dd.close();
      FileReader filetmp = new FileReader(pathtmp);
      input = new BufferedReader(filetmp);
      FileWriter escribedd = new FileWriter(path0, true);
      String linea;
      while ((linea = input.readLine()) != null) {
        String[] personas = linea.split("/");
        for (String persona : personas) {
          String[] creaAdanEva = persona.split(",");
            Persona sujetoprueba = new Persona();
            sujetoprueba.setNombre(creaAdanEva[0]);
            sujetoprueba.setApellido(creaAdanEva[1]);
            sujetoprueba.setNumCuenta(Long.parseLong(creaAdanEva[2].trim()));
            escribedd = new FileWriter(path0, true);
            escribedd.write(sujetoprueba.toString());
            escribedd.close();
          } // Fin foreach de cada sujeto
      } // fin while de lectura
      input.close();
      f.delete();
    } catch(Exception e) {
      e.printStackTrace();
    } // Fin del try-catch

  } // Fin del metodo

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
   * Para poder escribir en los archivos g1 o g2, alternando cada bloque que se tiene.
   * @param personas La lista de personas a escribir.
   * @param alterna El criterio para decidir si se escriben en g1 o g2.
   * @param path1 La ruta del archivo g1.
   * @param path2 La ruta del archivo g2.
   * @param path1itera La ruta del archivo g1_Iteraciones.
   * @param path2itera La ruta del archivo g2_Iteraciones.
   * @return Una lista de personas, que contiene a la ultima persona leida pero que es menor al ultimo del bloque.
   */
  public List<Persona> escribeArchivog1g2(List<Persona> personas,int alterna,String path1,String path2, String path1itera, String path2itera) {
    if (alterna % 2 != 0) {
      try { // Para archivo g1
        FileWriter g1 = new FileWriter(path1, true);
        BufferedWriter escribe = new BufferedWriter(g1);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        g1.close();
        FileWriter g1itera = new FileWriter(path1itera, true);
        escribe = new BufferedWriter(g1itera);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        g1itera.close();
        personas.clear();
        
      } catch (IOException e) {
        System.out.println("Algo anda mal");
      }
      return personas;
    } else {
      try { // Para archivo g2
        FileWriter g2 = new FileWriter(path2, true);
        BufferedWriter escribe = new BufferedWriter(g2);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        g2.close();
        FileWriter g2itera = new FileWriter(path2itera, true);
        escribe = new BufferedWriter(g2itera);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        g2itera.close();
        personas.clear();
        
      } catch (IOException e) {
        System.out.println("Algo anda mal");
      } // Fin try-catch
      return personas;
    } // Fin else que imprime en archivo g2
  } // Fin del metodo

  /**
   * Este metodo sirve para escribir en el archivo, al ultimo bloque leido desde dd.
   * @param personas La lista de personas a escribir en el archivo.
   * @param alterna El criterio para decidir si escribe en g1 o g2.
   * @param path1 La ruta del archivo g1.
   * @param path2 La ruta del archivo g2.
   * @param path1itera La ruta del archivo g1_Iteraciones.
   * @param path2itera La ruta del archivo g2_Iteraciones.
   */
  public void escribeArchivoFinal(List<Persona> personas,int alterna,String path1,String path2, String path1itera, String path2itera) {
    if (alterna % 2 != 0) {
      try { // Para archivo g1
        FileWriter g1 = new FileWriter(path1, true);
        BufferedWriter escribe = new BufferedWriter(g1);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        g1.close();
        FileWriter g1itera = new FileWriter(path1itera, true);
        escribe = new BufferedWriter(g1itera);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        g1itera.close();
        personas.clear();
      } catch (IOException e) {
        System.out.println("Algo anda mal");
      }
    } else {
      try { // Para archivo g2
        FileWriter g2 = new FileWriter(path2, true);
        BufferedWriter escribe = new BufferedWriter(g2);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        g2.close();
        FileWriter g2itera = new FileWriter(path2itera, true);
        escribe = new BufferedWriter(g2itera);
        for (Persona per1 : personas) {
          escribe.write(per1.impArchivoAux());
        }
        escribe.newLine();
        escribe.close();
        g2itera.close();
        personas.clear();
      } catch (IOException e) {
        System.out.println("Algo anda mal");
      } // Fin try-catch
    } // Fin else que imprime en archivo g2
  } // Fin del metodo

/**
   * Este metodo sirve para eliminar los archivos g1 y g2.
   * @param path1 La ruta del archivo g1.
   * @param path2 La ruta del archivo g2.
   */
  public void borrarG(String path1, String path2){
    File g1 = new File(path1);
    File g2 = new File(path2);
    g1.delete();
    g2.delete();

  }
  

}
