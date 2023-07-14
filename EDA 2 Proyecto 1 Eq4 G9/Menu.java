import java.util.Scanner;

/**	La clase Menu es una estructura de datos compuesta que representa a un Menu seleccionable.
	Sus atributos son:
	String texto - todo el menu con titulo y las opciones disponibles.
	int opciones - el numero de opciones que tiene el menu, incluyendo salir/regresar.

	Por medio de estos atributos las instancia de la clase menu pueden mostrar el menu y solicitar
	al usuario que ingrese la opcion deseada. Siempre dentro del rango permitido.
	@author David Calderon
	@author Humberto Hernandez
	@author Inaky Ordiales
*/
public class Menu{
	
	String texto;
	private int opciones;

	/**	El unico constructor de la clase Menu ya que no cuenta con el vacio.
		recibe el que sera el texto del menu como String textoMenu y el numero de opciones como
		int opciones. A traves de los metodos de acceso setters establece sus atributos. Si el setter
		para opciones devuelve falso, se despliega mensaje de error.
		@param textoMenu es todo el texto del menu que quieres que se despliegue en pantalla.
		@param opciones es el numero de opciones que se pueden elegir en el menu a desplegar.
	*/
	public Menu(String textoMenu, int opciones){
		setTexto(textoMenu);
		if (!setOpciones(opciones))
			System.out.println(" El numero de opciones probablemente sea incorrecto.");
	}
	/**	El metodo setter setTexto recibe como parametro una String que debe ser todo el texto del menu
		y lo asigna al atributo.
		@param texto es todo el texto del menu que quieres que se despliegue en pantalla.
	*/
	public void setTexto(String texto){
		this.texto = texto;
	}
	/** Este es un metodo getter de un atributo del objeto Menu.
		@return regresa el atributo texto del objeto Menu.
	*/
	public String getTexto(){
		return texto;
	}

	/**	El metodo setter setOpciones recibe como parametro un entero que debe ser el numero de opciones
		disponibles en el menu y regresa un booleano. 
		En caso contrario asigna correctamente el parametro entero al atributo opciones.
		@param opciones el numero de opciones que tienen el menu.
		@return Si las opciones son menores o iguales a uno, regresara false porque los menus de una opcion no hacen sentido.
	*/
	public boolean setOpciones(int opciones){
		if (opciones>1){
			this.opciones = opciones;
			return true;
		}
		return false;
	}
	/** Este es un metodo getter de un atributo del objeto Menu.
		@return regresa el atributo opciones del objeto Menu.
	*/
	public int getOpciones(){
		return opciones;
	}

	/**	El metodo mostrarMenu despliega en pantalla el texto del atributo texto del menu junto con un formato predefinido.
		Ademas al final imprimer 'OPCION: ' para que el usuario pueda elegir.
		Este metodo regresa un número entero correspondiente a la opcion que se eligio. Para hacer este return manda a llamar
		al metodo obtenerOpcion() y regresa el entero que a su vez obtenerOpcion regrese.
		@return el numero de opcion elegida.
	*/
	public int mostrarMenu(){
		System.out.println("\n\n\n"+getTexto());
		System.out.print("      OPCION: ");
		return obtenerOpcion();
	}
	/**	El metodo obtenerOpcion es el encargado de recibir y verificar el número de opción que fue elegido por el usuario.
		Dentro cuenta con una implementación de un ciclo que mientras el valor ingresado por el usuario no entre en el
		rango de 1 hasta el atributo opciones, manda mensaje del error a pantalla y llama recursivamente al metodo mostrarMenu()
		para empezar todo el proceso de nuevo. Nunca se va a regresar un entero fuera del rango de opciones.
		@return el numero de opcion del menu que el usuario ha elegido.
	*/
	public int obtenerOpcion(){
		Scanner sc = new Scanner(System.in);

		int opcion = sc.nextInt();
		while (opcion<1 || opcion>opciones){
			System.out.println("\n\t\t Opcion NO VALIDA.");
			return mostrarMenu();
		}
		return opcion;
	}

}