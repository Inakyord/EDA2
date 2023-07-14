package polifase;

/**
 * Clase principal para polifase.
 * @author David Calderon Jimenez
 * @author Humberto Ignacio Hernandez Olvera
 * @author I\u00f1aky Ordiales Caballero
 * @version 1.5
 * @since 29/Octubre/2020
 */

public class Polifase {
  public void ordenar(String path, boolean choice) {
    // Si es true es para apellido, false es para nombre
    boolean porApellido = choice;
    System.out.println("Ordenamiento por polifase con bloques de 4 elementos\n");
    PolifaseMetodos ordenamiento = new PolifaseMetodos();
    ordenamiento.borraArchivos(porApellido);
    // en 0 es path de archivo original, en 1 es de g1 y en 2 es de g2
    String[] pathsAux = new String[5];
    ordenamiento.creaArchivos(porApellido); 
    // Se copia el archivo a un nuevo directorio para mantener el original intacto
    //String path0 = ordenamiento.clonOriginal(path, porApellido); 

    pathsAux = ordenamiento.ordenaI1(path, porApellido); // {g0, g1, g2, g1ITERA, g2ITERA}
    int lineas = ordenamiento.comprueba(path);
    System.out.println("Se han detectado " + lineas + " registros a ordenar.\n");

    /*
    El caso para terminar el ciclo de ordenamiento, es cuando el archivo original (g0) contiene una sola línea
    */
    while (lineas > 1) {
      // Modifica g0. {g0, g1, g2, g1ITERA, g2ITERA}
      pathsAux = ordenamiento.ordenaI2(pathsAux[0] ,pathsAux[1], pathsAux[2], porApellido, pathsAux[3], pathsAux[4]);
       // Crea los bloques en g1 y g2 para proceder a ordenarlos despues
      pathsAux = ordenamiento.ordenaI3(pathsAux[0], pathsAux[1], pathsAux[2], pathsAux[3], pathsAux[4]);
      // Para revisar cuantas lineas tiene el archivo g0 
      lineas = ordenamiento.comprueba(pathsAux[0]);
      // Para revisar cuantas lineas tiene el archivo g1
      int lineas1 = ordenamiento.comprueba(pathsAux[1]);
      // Para revisar cuantas lineas tiene el archivo g2
      int lineas2 = ordenamiento.comprueba(pathsAux[2]);
      if (lineas1 == 1 && lineas2 == 1) {
        // Modifica g0. {g0, g1, g2}
        pathsAux = ordenamiento.ordenaI2(pathsAux[0] ,pathsAux[1], pathsAux[2], porApellido, pathsAux[3], pathsAux[4]); 
        break;
      }
    }
    // Se procede a regresar el archivo a su formato original.
    ordenamiento.regresaAFormato(pathsAux[0], porApellido);  
    ordenamiento.borrarG(pathsAux[1],pathsAux[2]);
    if (porApellido) {
      System.out.println("El archivo ordenado se encuentra en porApellido/g0.txt");
    } else {
      System.out.println("El archivo ordenado se encuentra en porNombre/g0.txt");
    }

  } // Fin del metodo
}
