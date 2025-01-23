class EmissionsCalculator:
    """
    A class to calculate emissions based on fuel consumption.
    """

    def __init__(self, fuel_consumption_tons, sulfur_content=0.5):
        """
        Initialize the EmissionsCalculator with fuel consumption and sulfur content.

        Parameters:
            fuel_consumption_tons (float): Fuel consumption in tons.
            sulfur_content (float): Sulfur content in the fuel (percentage, default is 0.5%).
        """
        self.fuel_consumption_tons = fuel_consumption_tons
        self.sulfur_content = sulfur_content

    def calculate_CO2_emissions(self):
        """
        Calculate CO2 emissions based on fuel consumption.

        Returns:
            float: CO2 emissions in tons.
        """
        EF_CO2 = 3.206  # CO2 emission factor (tons per ton)
        return self.fuel_consumption_tons * EF_CO2

    def calculate_SO2_emissions(self):
        """
        Calculate SO2 emissions based on fuel consumption and sulfur content.

        Returns:
            float: SO2 emissions in tons.
        """
        EF_SO2 = self.sulfur_content * 0.02  # SO2 emission factor depends on sulfur content (tons per ton)
        return self.fuel_consumption_tons * EF_SO2

    def calculate_NOx_emissions(self):
        """
        Calculate NOx emissions based on fuel consumption.

        Returns:
            float: NOx emissions in tons.
        """
        EF_NOx = 0.07  # NOx emission factor (tons per ton, average value)
        return self.fuel_consumption_tons * EF_NOx

if __name__ == "__main__":
    # Example usage
    fuel_consumption = 60  # Fuel consumption in tons
    sulfur_content = 0.5  # Sulfur content in percentage

    calculator = EmissionsCalculator(fuel_consumption, sulfur_content)

    CO2_emissions = calculator.calculate_CO2_emissions()
    SO2_emissions = calculator.calculate_SO2_emissions()
    NOx_emissions = calculator.calculate_NOx_emissions()

    print("Emissions based on fuel consumption:")
    print(f"CO2: {CO2_emissions:.2f} tons")
    print(f"SO2: {SO2_emissions:.2f} tons")
    print(f"NOx: {NOx_emissions:.2f} tons")