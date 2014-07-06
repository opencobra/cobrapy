{%- extends 'display_priority.tpl'-%}

{% block input %}
{{ cell.input }}
{% endblock input %}

{%- block pyout -%}
{%- block data_priority scoped -%}
{{ super() }}
{%- endblock -%}
{%- endblock pyout -%}

{% block data_text scoped -%}
# Output:
{{output.text | comment_lines }}
{% endblock data_text %}

{% block stream -%}
# Prints:
{{output.text | trim | comment_lines }}
{% endblock stream %}

{% block markdowncell scoped %}
{{ cell.source | wrap_text(77) | comment_lines }}
{% endblock markdowncell %}

{% block headingcell scoped %}
{{ '#' * cell.level }}{{ cell.source | replace('\n', ' ') | comment_lines }}
{% endblock headingcell %}
