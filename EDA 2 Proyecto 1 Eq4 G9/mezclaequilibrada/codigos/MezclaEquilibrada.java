// import generador.*;
// import java.io.FileWriter;
// import java.io.IOException;
// import java.io.File;
// import java.io.FileNotFoundException;
// import java.util.Scanner;

package mezclaequilibrada.codigos;

/**
 * Clase principal para mezcla equilibrada.
 * @author David Calderon Jimenez
 * @author Humberto Ignacio Hernandez Olvera
 * @author I\u00f1aky Ordiales Caballero
 * @version 1.9
 * @since 22/Octubre/2020
 */

public class MezclaEquilibrada {
  public void ordenar(String path, boolean choice) {
    // Si es true es para apellido, false es para nombre
    boolean porApellido = choice;
    // System.out.println("Se creara una copia del archivo original, para no modificarlo y asi usarlo en mas pruebas.\n");

    MezclaEquilibradaMetodos ordenamiento = new MezclaEquilibradaMetodos();
    ordenamiento.limpiaArchivos(porApellido);
    // en 0 es path de archivo original, en 1 es de F1 y en 2 es de F2
    String[] pathsAux = new String[5];
    ordenamiento.creaArchivos(porApellido);
    if (porApellido) {
      ordenamiento.llamadaAp++;
    } else {
      ordenamiento.llamadaNom++;
    }
    // Se copia el archivo a un nuevo directorio para mantener el original intacto
    // String path0 = ordenamiento.clonOriginal(path, porApellido);
    String path0 = path;

    pathsAux = ordenamiento.ordenaI1(path0, porApellido); // {archivo a modificar, F1, F2, F1ITERA, F2ITERA}
    int lineas = ordenamiento.comprueba(path0);
    System.out.println("\n\t En el archivo "+path0+" se encontraron "+lineas+" registros.\n");

    /*
    El caso para terminar el ciclo de ordenamiento, es cuando el archivo original (F0) contiene una sola línea
    */
    while (lineas > 1) {
      // {archivo a modificar, F1, F2, F1ITERA, F2ITERA}
      pathsAux = ordenamiento.ordenaI2(pathsAux[0] ,pathsAux[1], pathsAux[2], porApellido, pathsAux[3], pathsAux[4]);
       // Crea los bloques en F1 y F2 para proceder a ordenarlos despues
      pathsAux = ordenamiento.ordenaI3(pathsAux[0], pathsAux[1], pathsAux[2], pathsAux[3], pathsAux[4]);
      // Para revisar cuantas lineas tiene el archivo a modificar
      lineas = ordenamiento.comprueba(pathsAux[0]);
      // Para revisar cuantas lineas tiene el archivo F1
      int lineas1 = ordenamiento.comprueba(pathsAux[1]);
      // Para revisar cuantas lineas tiene el archivo F2
      int lineas2 = ordenamiento.comprueba(pathsAux[2]);
      if (lineas1 == 1 && lineas2 == 1) {
        // {archivo a modificar, F1, F2, F1ITERA, F2ITERA}
        pathsAux = ordenamiento.ordenaI2(pathsAux[0] ,pathsAux[1], pathsAux[2], porApellido, pathsAux[3], pathsAux[4]);
        break;
      }
    }
    // Se procede a regresar el archivo a su formato original.
    ordenamiento.regresaAFormato(pathsAux[0], porApellido);
    if (porApellido) {
      System.out.println("\tLos archivos con las iteraciones se encuentran en: ");
      System.out.println("\t\t"+pathsAux[3]);
      System.out.println("\t\t"+pathsAux[4]);
    } else {
      System.out.println("\tLos archivos con las iteraciones se encuentran en: ");
      System.out.println("\t\t"+pathsAux[3]);
      System.out.println("\t\t"+pathsAux[4]);
    }
    ordenamiento.borraArchivos(pathsAux[1], pathsAux[2]);

  } // Fin del metodo
}
