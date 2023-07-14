package polifase;

/**
 * Esta clase sirve para crear a la persona que se lee del archivo.
 * @author David Calderon Jimenez
 * @author Humberto Ignacio Hernandez Olvera
 * @author I\u00f1aky Ordiales Caballero
 * @version 1.0
 * @since 22/Octubre/2020
 */
public class Persona {
  String nombre;
  String apellido;
  long num_cuenta;

  /**
   * Constructor vacio para persona.
   */
  public Persona() {

  }

  /**
   * Constructor personalizado para persona.
   * @param nombre El nombre o nombres de la persona.
   * @param apellido Los apellidos de la persona.
   * @param num_cuenta El numero de cuenta de la persona.
   */
  public Persona(String nombre, String apellido, long num_cuenta) {
    this.nombre = nombre;
    this.apellido = apellido;
    this.num_cuenta = num_cuenta;
  }

  /**
   * Metodo setter para el nombre.
   * @param nombre El nombre o nombres de la persona.
   */
  public void setNombre(String nombre){
    this.nombre = nombre;
  }

  /**
   * Metodo setter para el apellido.
   * @param apellido Los apellidos de la persona.
   */
  public void setApellido(String apellido){
    this.apellido = apellido;
  }

  /**
   * Metodo setter para el numero de cuenta.
   * @param num_cuenta El numero de cuenta de la persona.
   */
  public void setNumCuenta(long num_cuenta) {
    this.num_cuenta = num_cuenta;
  }

  /**
   * Metodo getter para el nombre.
   * @return El nombre o nombres de la persona.
   */
  public String getNombre() {
    return this.nombre;
  }

  /**
   * Metodo getter para el apellido.
   * @return Los apellidos de la persona.
   */
  public String getApellido() {
    return this.apellido;
  }

  /**
   * Metodo getter para el numero de cuenta.
   * @return El numero de cuenta de la persona.
   */
  public long getNumCuenta() {
    return this.num_cuenta;
  }

  /**
   * Metodo que estandariza la escritura de los bloques en los archivos durante el ordenamiento.
   * @return La cadena que contiene a todas las personas que representaran una linea en el archivo.
   */
  public String impArchivoAux() {
    return this.nombre + "," + this.apellido + ","  + this.num_cuenta + "/";
  }

  /**
   * Sobreescritura del metodo toString(), que ademas sirve para imprimir en el formato original del archivo.
   */
  @Override
  public String toString() {
    return this.nombre + "," + this.apellido + ","  + this.num_cuenta + "\n";
  }

}
